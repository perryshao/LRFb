#!/usr/bin/env python3
"""Compare a formatting-only revision with a pre-polish LRFb tar.gz snapshot.

Requires miss-hit==0.9.44 and clang++. This checks MATLAB syntax trees and
significant tokens, C/C++ raw tokens, and byte identity of vendored code and
binary assets. It does not execute MATLAB or compile/link the MEX libraries.
"""

import argparse
import hashlib
import io
import json
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

from miss_hit_core.config import Config
from miss_hit_core.errors import Message_Handler
from miss_hit_core.m_lexer import MATLAB_Lexer, Token_Buffer
from miss_hit_core.m_parser import MATLAB_Parser


def matlab_signature(content, filename):
    """Keep expression structure and code tokens, ignoring comment layout."""
    content = content.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    config = Config()
    config.style_rules = set()
    messages = Message_Handler("style")
    messages.register_file(filename)
    messages.show_style = False
    lexer = MATLAB_Lexer(config.language, messages, content, filename)
    buffer = Token_Buffer(lexer, config)
    tree = MATLAB_Parser(messages, buffer, config).parse_file()
    output = io.StringIO()
    tree.pp_node(output)
    tokens = []
    for token in buffer.tokens:
        if token.anonymous or token.kind in {"COMMENT", "CONTINUATION", "NEWLINE"}:
            continue
        if token.kind == "COMMA" and token.fix.statement_terminator:
            # Replacing a statement-ending comma with a newline preserves output.
            continue
        tokens.append((token.kind, token.value))
    return output.getvalue(), tokens


def cpp_signature(content, filename, clang, directory):
    """Lex without preprocessing, so MATLAB headers are not required."""
    path = directory / Path(filename).name
    path.write_bytes(content)
    result = subprocess.run(
        [clang, "-x", "c++", "-fsyntax-only", "-Xclang", "-dump-raw-tokens", str(path)],
        capture_output=True,
        check=True,
    )
    dump = result.stderr.decode("utf-8", errors="replace")
    pattern = r"^(\w+) '(.*?)'\s*(?:\[[^\]]*\]\s*)*Loc=<[^\n]*>$"
    tokens = []
    for match in re.finditer(pattern, dump, re.MULTILINE | re.DOTALL):
        kind, spelling = match.groups()
        if kind == "comment":
            continue
        if kind == "unknown" and not spelling.strip():
            continue
        tokens.append((kind, spelling))
    if not tokens:
        raise RuntimeError(f"No Clang tokens parsed for {filename}")
    return tokens


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--json", type=Path, help="Write the verification results")
    args = parser.parse_args()
    clang = shutil.which("clang++")
    if not clang:
        parser.error("clang++ is required for C/C++ token comparison")

    results = []
    with tarfile.open(args.snapshot) as archive, tempfile.TemporaryDirectory() as temp:
        scratch = Path(temp)
        for member in archive.getmembers():
            if not member.isfile():
                continue
            parts = Path(member.name).parts
            if not parts or parts[0] != "LRFb" or ".." in parts:
                continue
            if any(part.startswith("._") for part in parts):
                # macOS tar synthesizes AppleDouble metadata entries from xattrs.
                continue
            relative = Path(*parts[1:])
            source = args.root / relative
            vendor = "thirdparty" in relative.parts
            native = relative.suffix in {".cpp", ".h", ".c"}
            matlab = relative.suffix == ".m"
            binary = relative.suffix.startswith(".mex")
            if not (vendor or native or matlab or binary):
                continue
            old = archive.extractfile(member).read()
            if not source.is_file():
                results.append({"file": str(relative), "check": "exists", "passed": False})
                continue
            new = source.read_bytes()
            if vendor or binary:
                check, passed = "byte_identity", old == new
            elif matlab:
                check = "matlab_ast_and_tokens"
                passed = matlab_signature(old, str(relative)) == matlab_signature(new, str(relative))
            else:
                check = "cpp_raw_tokens"
                passed = cpp_signature(old, str(relative), clang, scratch) == cpp_signature(
                    new, str(relative), clang, scratch
                )
            results.append({
                "file": str(relative),
                "check": check,
                "passed": passed,
                "before_sha256": hashlib.sha256(old).hexdigest(),
                "after_sha256": hashlib.sha256(new).hexdigest(),
            })
    if not results:
        raise SystemExit("Snapshot contains no matching LRFb source or binary files")
    failures = [row for row in results if not row["passed"]]
    if args.json:
        args.json.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"Verified {len(results)} files; {len(failures)} differences requiring review.")
    for row in failures:
        print(f"FAIL {row['check']}: {row['file']}")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
