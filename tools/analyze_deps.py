#!/usr/bin/env python3
"""
analyze_deps.py -- build a Makefile + source-code dependency graph as JSON.

Generic two-level analyzer:
  1) Make level: parses one or more Makefiles (GNU Make style: explicit
     rules "target: prereqs" + tab-indented recipes, simple "VAR = value"
     variables with $(VAR) substitution). Detects:
       - which targets are final executables (link recipes: contain "-o $@"
         and no "-c" flag)
       - which .o targets map to which source file (recipes containing "-c
         <source>")
       - recursive make invocations in recipes ($(MAKE)/make -C/-f) as
         edges to other Makefiles
  2) Source level: parses each referenced source file (Fortran free-form,
     case-insensitive) for MODULE/PROGRAM/SUBROUTINE/FUNCTION definitions,
     USE statements, and CALL statements + bare-function-call heuristics,
     to build a routine-level call graph.

This is NOT a full Make or Fortran parser -- no conditionals, no Make
functions ($(wildcard ...) etc.), no pattern rules (%.o: %.f90), no
Fortran preprocessor. It targets the common style of explicit-rule
Makefiles + free-form Fortran used in HPC/scientific codebases. Safe to
re-run any time the source changes; purely read-only (never touches the
inspected repository).

Usage (run from the repository root):
    python3 tools/analyze_deps.py OUTPUT.json Makefile [Makefile2 ...]

Each Makefile argument is parsed independently (and any Makefile it
recursively invokes is followed automatically); pass several to also
cover related-but-separate builds (e.g. a legacy/Makefile that is never
called BY the main Makefile but shares source files with it). For this
repository:

    python3 tools/analyze_deps.py depgraph.json Makefile legacy/Makefile
"""
import json
import os
import re
import sys
from collections import defaultdict

# --------------------------------------------------------------------------
# Part 1: Makefile parsing
# --------------------------------------------------------------------------

VAR_ASSIGN_RE = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)\s*(:=|\+=|=)\s*(.*)$')
RULE_RE = re.compile(r'^([^\t#:][^:]*):(?!=)\s*(.*)$')
VARREF_RE = re.compile(r'\$\(([A-Za-z_][A-Za-z0-9_]*)\)')
RECURSIVE_MAKE_RE = re.compile(
    r'\$\(MAKE\)|(?<![A-Za-z_])make\b', re.IGNORECASE)
DASH_C_RE = re.compile(r'(?:^|\s)-C\s+(\S+)')
DASH_F_RE = re.compile(r'(?:^|\s)-f\s+(\S+)')
SOURCE_EXTS = ('.f90', '.f', '.F90', '.F', '.c', '.cc', '.cpp', '.cxx')


def read_logical_lines(path):
    """Read a Makefile, joining backslash line-continuations, stripping
    comments (naively: a '#' anywhere outside a recipe line's shell
    quoting -- good enough for these Makefiles), and keeping track of
    which lines are recipe lines (start with a literal TAB)."""
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        raw = f.read()
    raw = raw.replace('\r\n', '\n')
    physical = raw.split('\n')

    logical = []  # list of (is_recipe, text)
    buf = ''
    buf_is_recipe = False
    for line in physical:
        is_recipe = line.startswith('\t')
        content = line[1:] if is_recipe else line
        if buf:
            content = buf + content
            buf = ''
        if content.rstrip().endswith('\\') and not is_recipe:
            buf = content.rstrip()[:-1]
            buf_is_recipe = is_recipe
            continue
        if not is_recipe:
            # strip comments outside recipes only (recipes are shell code,
            # '#' there is usually a shell comment too, but keep them
            # verbatim since we only regex-search them, never execute)
            hashpos = content.find('#')
            if hashpos >= 0:
                content = content[:hashpos]
        logical.append((is_recipe, content))
    if buf:
        logical.append((buf_is_recipe, buf))
    return logical


def expand_vars(text, variables, _depth=0):
    if _depth > 20 or '$(' not in text:
        return text
    def repl(m):
        name = m.group(1)
        return variables.get(name, '')
    new_text = VARREF_RE.sub(repl, text)
    if new_text == text:
        return text
    return expand_vars(new_text, variables, _depth + 1)


def parse_makefile(path):
    """Returns dict: variables, rules {target: {'prereqs':[...],
    'recipe':[...]}}, recursive_make_calls (list of resolved paths/hints)."""
    lines = read_logical_lines(path)
    variables = {}
    rules = {}
    order = []  # target definition order, for stable UI ordering
    recursive_calls = []

    current_targets = None
    for is_recipe, content in lines:
        if is_recipe:
            if current_targets:
                expanded = expand_vars(content, variables)
                for t in current_targets:
                    rules[t]['recipe'].append(expanded)
                if RECURSIVE_MAKE_RE.search(expanded):
                    hint = expanded.strip()
                    cdir = DASH_C_RE.search(hint)
                    ffile = DASH_F_RE.search(hint)
                    recursive_calls.append({
                        'from_recipe': hint,
                        'dir': cdir.group(1) if cdir else None,
                        'file': ffile.group(1) if ffile else None,
                    })
            continue

        stripped = content.strip()
        if not stripped:
            current_targets = None
            continue

        vm = VAR_ASSIGN_RE.match(stripped)
        if vm and ':' not in stripped.split('=')[0]:
            name, op, val = vm.groups()
            val = expand_vars(val.strip(), variables)
            if op == '+=' and name in variables:
                variables[name] = variables[name] + ' ' + val
            else:
                variables[name] = val
            current_targets = None
            continue

        rm = RULE_RE.match(stripped)
        if rm:
            targets_raw, prereqs_raw = rm.groups()
            targets_raw = expand_vars(targets_raw, variables)
            prereqs_raw = expand_vars(prereqs_raw, variables)
            targets = targets_raw.split()
            prereqs = prereqs_raw.split()
            current_targets = targets
            for t in targets:
                if t not in rules:
                    rules[t] = {'prereqs': [], 'recipe': []}
                    order.append(t)
                for p in prereqs:
                    if p not in rules[t]['prereqs']:
                        rules[t]['prereqs'].append(p)
            continue
        current_targets = None

    return {
        'path': path,
        'variables': variables,
        'rules': rules,
        'order': order,
        'recursive_calls': recursive_calls,
    }


def is_link_recipe(recipe_lines):
    text = ' '.join(recipe_lines)
    return ('-o $@' in text or '-o  $@' in text) and re.search(r'(?<!-c)\s-c\s', text) is None \
        and not re.search(r'\s-c\s', text)


def find_source_in_recipe(recipe_lines):
    text = ' '.join(recipe_lines)
    for tok in text.split():
        tok_clean = tok.strip()
        if tok_clean.endswith(SOURCE_EXTS):
            return tok_clean
    return None


def classify_makefile(mk):
    """Decide which targets are 'executables' (link rules) vs 'objects'
    (compile rules) vs 'aggregate/phony' (neither -- pure grouping
    targets like 'all'/'examples'/'clean')."""
    executables = []
    objects = {}   # target -> source file (or None if not determined)
    aggregates = []

    for t in mk['order']:
        rule = mk['rules'][t]
        recipe = rule['recipe']
        if not recipe:
            if rule['prereqs']:
                aggregates.append(t)
            continue
        if t.endswith('.o'):
            src = find_source_in_recipe(recipe)
            objects[t] = src
            continue
        if is_link_recipe(recipe):
            executables.append(t)
        else:
            # has a recipe but isn't an obvious link/compile step
            # (e.g. a 'clean' rule) -- ignore for the dependency graph
            pass
    return executables, objects, aggregates


def resolve_prereq_sources(target, mk, objects, base_dir, visited=None):
    """Flatten a target's prerequisite tree down to the set of actual
    source files (walks through intermediate .o/aggregate targets).
    Paths are re-rooted to be relative to the analyzer's own working
    directory, since a Makefile's own paths are relative to ITS
    directory (e.g. legacy/Makefile's "src/foo.f90" means
    "legacy/src/foo.f90" seen from the repo root)."""
    if visited is None:
        visited = set()
    if target in visited:
        return set()
    visited.add(target)
    sources = set()

    def reroot(p):
        return os.path.normpath(os.path.join(base_dir, p)) if base_dir not in ('', '.') else os.path.normpath(p)

    rule = mk['rules'].get(target)
    if rule is None:
        if target.endswith(SOURCE_EXTS):
            sources.add(reroot(target))
        return sources
    if target in objects and objects[target]:
        sources.add(reroot(objects[target]))
    for p in rule['prereqs']:
        if p.endswith(SOURCE_EXTS):
            sources.add(reroot(p))
        else:
            sources |= resolve_prereq_sources(p, mk, objects, base_dir, visited)
    return sources


def gather_makefiles(seed_paths):
    """Parse every seed Makefile; follow recursive $(MAKE)/-C/-f hints to
    discover related Makefiles automatically, without ever executing
    anything."""
    seen = {}
    queue = [os.path.normpath(p) for p in seed_paths]
    while queue:
        p = queue.pop(0)
        if p in seen or not os.path.isfile(p):
            continue
        mk = parse_makefile(p)
        seen[p] = mk
        base_dir = os.path.dirname(p) or '.'
        for call in mk['recursive_calls']:
            candidate_dir = call['dir']
            candidate_file = call['file']
            if candidate_file:
                cand = candidate_file if os.path.isabs(candidate_file) \
                    else os.path.normpath(os.path.join(
                        base_dir, candidate_dir or '', candidate_file))
                if os.path.isfile(cand):
                    queue.append(cand)
                    continue
            if candidate_dir:
                cand = os.path.normpath(os.path.join(base_dir, candidate_dir, 'Makefile'))
                if os.path.isfile(cand):
                    queue.append(cand)
    return seen


# --------------------------------------------------------------------------
# Part 2: Fortran source parsing
# --------------------------------------------------------------------------

MODULE_RE = re.compile(r'^\s*MODULE\s+([A-Za-z_][A-Za-z0-9_]*)\s*$', re.IGNORECASE)
END_MODULE_RE = re.compile(r'^\s*END\s*MODULE\b', re.IGNORECASE)
PROGRAM_RE = re.compile(r'^\s*PROGRAM\s+([A-Za-z_][A-Za-z0-9_]*)\s*$', re.IGNORECASE)
END_PROGRAM_RE = re.compile(r'^\s*END\s*PROGRAM\b', re.IGNORECASE)
SUBROUTINE_RE = re.compile(
    r'^\s*(?:RECURSIVE\s+|PURE\s+|ELEMENTAL\s+)*SUBROUTINE\s+([A-Za-z_][A-Za-z0-9_]*)', re.IGNORECASE)
END_SUBROUTINE_RE = re.compile(r'^\s*END\s*SUBROUTINE\b', re.IGNORECASE)
FUNCTION_RE = re.compile(
    r'^\s*(?!END\b)(?:[A-Za-z_][A-Za-z0-9_()*]*\s+)*?(?:RECURSIVE\s+|PURE\s+|ELEMENTAL\s+)*FUNCTION\s+([A-Za-z_][A-Za-z0-9_]*)',
    re.IGNORECASE)
END_FUNCTION_RE = re.compile(r'^\s*END\s*FUNCTION\b', re.IGNORECASE)
END_GENERIC_RE = re.compile(r'^\s*END\s*$', re.IGNORECASE)
CONTAINS_RE = re.compile(r'^\s*CONTAINS\s*$', re.IGNORECASE)
USE_RE = re.compile(r'^\s*USE\b\s*(?:,\s*INTRINSIC)?\s*(?:::)?\s*([A-Za-z_][A-Za-z0-9_]*)', re.IGNORECASE)
CALL_RE = re.compile(r'\bCALL\s+([A-Za-z_][A-Za-z0-9_]*)', re.IGNORECASE)
IDENT_CALL_RE = re.compile(r'\b([A-Za-z_][A-Za-z0-9_]*)\s*\(')
INTERFACE_RE = re.compile(r'^\s*INTERFACE\b', re.IGNORECASE)
END_INTERFACE_RE = re.compile(r'^\s*END\s*INTERFACE\b', re.IGNORECASE)

FORTRAN_KEYWORDS_NOT_CALLS = set(w.upper() for w in (
    'if', 'then', 'else', 'elseif', 'endif', 'do', 'enddo', 'end', 'while',
    'select', 'case', 'default', 'return', 'stop', 'continue', 'format',
    'real', 'integer', 'complex', 'logical', 'character', 'type', 'kind',
    'intent', 'in', 'out', 'inout', 'allocatable', 'pointer', 'parameter',
    'dimension', 'allocate', 'deallocate', 'write', 'read', 'open', 'close',
    'implicit', 'none', 'contains', 'module', 'program', 'subroutine',
    'function', 'use', 'only', 'result', 'recursive', 'pure', 'elemental',
    'public', 'private', 'interface', 'abstract', 'procedure', 'external',
    'optional', 'target', 'save', 'common', 'equivalence', 'data',
    'present', 'associated', 'allocated', 'size', 'shape', 'min', 'max',
    'abs', 'sqrt', 'real', 'aimag', 'cmplx', 'int', 'nint', 'mod', 'sum',
    'merge', 'trim', 'adjustl', 'adjustr', 'len', 'len_trim', 'reshape',
    'transpose', 'matmul', 'huge', 'tiny', 'kind', 'selected_real_kind',
))


def strip_fortran_comment(line):
    # crude: '!' outside of a quoted string starts a comment
    in_s = None
    for i, ch in enumerate(line):
        if in_s:
            if ch == in_s:
                in_s = None
        elif ch in ('"', "'"):
            in_s = ch
        elif ch == '!':
            return line[:i]
    return line


def join_continuations(lines):
    out = []
    buf = ''
    for line in lines:
        code = strip_fortran_comment(line)
        rs = code.rstrip()
        if rs.endswith('&'):
            buf += rs[:-1]
            continue
        if buf:
            out.append(buf + code)
            buf = ''
        else:
            out.append(code)
    if buf:
        out.append(buf)
    return out


def parse_fortran_file(path, all_known_names):
    """Returns dict: modules -> [names], routines -> [ {name, kind,
    module_or_program, line, calls:set, uses:set} ]."""
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            raw_lines = f.readlines()
    except OSError:
        return {'modules': [], 'routines': []}

    lines = join_continuations([l.rstrip('\n') for l in raw_lines])

    modules = []
    routines = []
    stack = []            # stack of ('module'|'program', name)
    current_routine = None
    interface_depth = 0
    module_uses = defaultdict(set)   # module/program name -> host-scope USEs,
                                      # inherited by every routine it CONTAINS

    def flush_routine():
        nonlocal current_routine
        if current_routine is not None:
            routines.append(current_routine)
            current_routine = None

    for lineno, line in enumerate(lines, start=1):
        s = line.strip()
        if not s:
            continue

        if INTERFACE_RE.match(s):
            interface_depth += 1
            continue
        if END_INTERFACE_RE.match(s):
            interface_depth = max(0, interface_depth - 1)
            continue
        if interface_depth > 0:
            continue  # skip interface bodies: declarations, not real calls

        m = MODULE_RE.match(s)
        if m:
            name = m.group(1)
            modules.append(name)
            stack.append(('module', name))
            continue
        if END_MODULE_RE.match(s):
            flush_routine()
            if stack and stack[-1][0] == 'module':
                stack.pop()
            continue

        m = PROGRAM_RE.match(s)
        if m:
            name = m.group(1)
            stack.append(('program', name))
            continue
        if END_PROGRAM_RE.match(s):
            flush_routine()
            if stack and stack[-1][0] == 'program':
                stack.pop()
            continue

        if CONTAINS_RE.match(s):
            continue

        m = SUBROUTINE_RE.match(s)
        if m:
            flush_routine()
            owner = stack[-1][1] if stack else None
            owner_kind = stack[-1][0] if stack else None
            current_routine = {
                'name': m.group(1).upper(), 'kind': 'subroutine',
                'owner': owner, 'owner_kind': owner_kind,
                'line': lineno, 'calls': set(), 'uses': set(),
            }
            stack.append(('subroutine', m.group(1).upper()))
            continue
        if END_SUBROUTINE_RE.match(s):
            if stack and stack[-1][0] == 'subroutine':
                stack.pop()
            flush_routine()
            continue

        m = FUNCTION_RE.match(s)
        if m:
            flush_routine()
            owner = stack[-1][1] if stack else None
            owner_kind = stack[-1][0] if stack else None
            current_routine = {
                'name': m.group(1).upper(), 'kind': 'function',
                'owner': owner, 'owner_kind': owner_kind,
                'line': lineno, 'calls': set(), 'uses': set(),
            }
            stack.append(('function', m.group(1).upper()))
            continue
        if END_FUNCTION_RE.match(s):
            if stack and stack[-1][0] == 'function':
                stack.pop()
            flush_routine()
            continue

        if END_GENERIC_RE.match(s) and stack and stack[-1][0] in ('subroutine', 'function'):
            stack.pop()
            flush_routine()
            continue

        um = USE_RE.match(s)
        if um:
            if current_routine is not None:
                current_routine['uses'].add(um.group(1).upper())
            elif stack:
                module_uses[stack[-1][1]].add(um.group(1).upper())
            continue

        if current_routine is None:
            continue

        for cm in CALL_RE.finditer(s):
            current_routine['calls'].add(cm.group(1).upper())
        for im in IDENT_CALL_RE.finditer(s):
            name = im.group(1).upper()
            if name in FORTRAN_KEYWORDS_NOT_CALLS:
                continue
            if name in all_known_names:
                current_routine['calls'].add(name)

    flush_routine()
    for r in routines:
        if r['owner'] in module_uses:
            r['uses'] |= module_uses[r['owner']]
        r['calls'] = sorted(r['calls'])
        r['uses'] = sorted(r['uses'])
    return {'modules': modules, 'routines': routines}


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    out_path = sys.argv[1]
    seed_makefiles = sys.argv[2:]

    makefiles = gather_makefiles(seed_makefiles)

    exec_entries = []
    all_sources = set()
    mk_summaries = []
    for path, mk in makefiles.items():
        executables, objects, aggregates = classify_makefile(mk)
        mk_summaries.append({
            'path': path,
            'n_targets': len(mk['order']),
            'n_executables': len(executables),
            'recursive_calls': mk['recursive_calls'],
        })
        base_dir = os.path.dirname(path)
        for exe in executables:
            srcs = sorted(resolve_prereq_sources(exe, mk, objects, base_dir))
            all_sources |= set(srcs)
            exec_entries.append({
                'name': exe, 'makefile': path, 'sources': srcs,
            })

    # First pass: collect every routine/module name across all sources so
    # the bare-identifier-call heuristic (IDENT_CALL_RE) has something to
    # match against.
    prelim = {}
    for src in sorted(all_sources):
        if not os.path.isfile(src):
            continue
        prelim[src] = parse_fortran_file(src, set())
    known_names = set()
    for src, data in prelim.items():
        for r in data['routines']:
            known_names.add(r['name'])

    files_out = {}
    routines_index = defaultdict(list)
    for src in sorted(all_sources):
        if not os.path.isfile(src):
            files_out[src] = {'modules': [], 'routines': [], 'missing': True}
            continue
        data = parse_fortran_file(src, known_names)
        files_out[src] = {'modules': data['modules'], 'routines': data['routines'], 'missing': False}
        for r in data['routines']:
            routines_index[r['name']].append({
                'file': src, 'line': r['line'], 'kind': r['kind'],
                'owner': r['owner'], 'owner_kind': r['owner_kind'],
            })

    # reverse call edges (callers) for convenience in the UI
    called_by = defaultdict(set)
    for src, data in files_out.items():
        for r in data['routines']:
            for callee in r['calls']:
                called_by[callee].add(f"{src}::{r['name']}")

    for src, data in files_out.items():
        for r in data['routines']:
            key = f"{src}::{r['name']}"
            r['called_by'] = sorted(called_by.get(r['name'], set()) - {key})

    result = {
        'makefiles': mk_summaries,
        'executables': sorted(exec_entries, key=lambda e: e['name']),
        'files': files_out,
        'routines_index': dict(routines_index),
    }
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(result, f, separators=(',', ':'), sort_keys=True)
    print(f"Wrote {out_path}: {len(exec_entries)} executables, "
          f"{len(all_sources)} source files, {len(known_names)} routines.")


if __name__ == '__main__':
    main()
