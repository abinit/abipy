#!/usr/bin/env python

import subprocess
import sys
import os

def get_changed_files():
    result = subprocess.run(
        #["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"],
        ["git", "diff", "--name-only", "--diff-filter=ACMR"],
        stdout=subprocess.PIPE, text=True
    )
    lines = [l.strip() for l in result.stdout.splitlines()]
    return [l for l in lines if l.endswith(".py")]

def get_test_files(files):
    # Adjust this logic to match your test naming scheme
    test_paths = []
    for path in files:
        dirname = os.path.dirname(path)
        basename = os.path.basename(path)
        test_path = os.path.join(dirname, "tests", f"test_{basename}")
        if os.path.exists(test_path):
            print(test_path)
            test_paths.append(test_path)

    return test_paths


def main():
    changed_files = get_changed_files()
    if len(changed_files) == 0:
        print("No staged Python files detected.")
        return 0

    test_files = get_test_files(changed_files)
    if not test_files:
        print("No modified test files found. Skipping test execution.")
        return 0

    print("Running pytest on modified test files:")
    for f in test_files:
        print(f" - {f}")

    cmd = "pytest -v " + " ".join(test_files)
    print("Running ", cmd)
    return os.system(cmd)

    #result = subprocess.run(["pytest -v", *test_files])
    #return result.returncode
    #return 0


if __name__ == "__main__":
    sys.exit(main())
