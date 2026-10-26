"""Validate JSON, then ask PARI/GP for an exact certificate verdict."""
import argparse
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

RATIONAL = re.compile(r"-?[0-9]+(?:/[1-9][0-9]*)?\Z", re.ASCII)
DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z", re.ASCII)


def object_without_duplicates(pairs):
    out = {}
    for key, value in pairs:
        if key in out:
            raise ValueError("duplicate JSON key")
        out[key] = value
    return out


def rational(value):
    if not isinstance(value, str) or RATIONAL.fullmatch(value) is None:
        raise ValueError("exact integer or fraction string required")
    return str(Fraction(value))


def decimal(value):
    if not isinstance(value, str) or DECIMAL.fullmatch(value) is None:
        raise ValueError("decimal claim required")
    result = Fraction(value)
    if result <= 0:
        raise ValueError("claim must be positive")
    return str(result)


def load(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=object_without_duplicates)
    if not isinstance(data, dict) or set(data) != {"n", "K", "tau0", "E", "c", "Gam", "claims"}:
        raise ValueError("invalid certificate keys")
    if type(data["n"]) is not int or data["n"] < 2 or type(data["K"]) is not int or data["K"] < 1:
        raise ValueError("invalid dimension or degree")
    values = [str(data["n"]), str(data["K"]), rational(data["tau0"])]
    for key in ("E", "c", "Gam"):
        if type(data[key]) is not list:
            raise ValueError("array required")
        values.append("[" + ",".join(rational(x) for x in data[key]) + "]")
    claims = data["claims"]
    if type(claims) is not dict or set(claims) != {"density", "hermite"}:
        raise ValueError("invalid claim keys")
    values += [decimal(claims["density"]), decimal(claims["hermite"])]
    return data, "[" + ",".join(values) + "]"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate")
    parser.add_argument("--gp", default=os.environ.get("PARI_GP", shutil.which("gp")))
    parser.add_argument("--timeout", type=float, default=3600)
    args = parser.parse_args()
    if not args.gp or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("supply --gp or PARI_GP; timeout must be positive")
    try:
        data, value = load(args.certificate)
        # The source path is generated locally; all certificate tokens were parsed above.
        source = (Path(__file__).parent / "verify.gp").resolve().as_posix()
        program = ('default(parisizemax,4000000000);\n'
                   + 'read(' + json.dumps(source) + ');\n'
                   + 'iferr(cc_verify(' + value + ');print("CERTIFICATE_ACCEPTED"),err,print("REJECT ",err);quit(1));\nquit(0);\n')
        proc = subprocess.run([args.gp, "-q", "-f"], input=program, text=True,
                              capture_output=True, timeout=args.timeout, cwd=Path(__file__).parent)
        if proc.returncode != 0 or proc.stdout.splitlines().count("CERTIFICATE_ACCEPTED") != 1:
            detail = (proc.stdout + proc.stderr).strip()
            raise ValueError("PARI rejection: " + detail[:1600])
        # GP reports parse errors on stderr without reliably failing its process.
        diagnostics = re.sub(r"(?m)^\s*\*\*\*\s+Warning: new maximum stack size = [0-9]+ \([0-9.]+ Mbytes\)\.\s*$", "", proc.stderr)
        diagnostics = re.sub(r"(?m)^\s*\*\*\*\s+polsturm: Warning: increasing stack size to [0-9]+\.\s*$", "", diagnostics)
        if "***" in diagnostics or "***" in proc.stdout:
            raise ValueError("PARI diagnostic indicates an incomplete computation: " + proc.stderr[:1800] + proc.stdout[:300])
        print(f'PASS n={data["n"]} K={data["K"]} density={data["claims"]["density"]} hermite={data["claims"]["hermite"]}')
        return 0
    except (OSError, ValueError, TypeError, OverflowError, subprocess.SubprocessError) as exc:
        print("REJECT " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
