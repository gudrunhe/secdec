if __name__ == "__main__":
    import os
    import sys
    from . import dirname
    argv = sys.argv[1:]
    if argv == ["--dirname"]:
        print(dirname)
    elif len(argv) >= 1:
        filename = os.path.join(dirname, "bin", argv[0])
        if os.path.isfile(filename):
            os.execv(filename, [filename] + argv[1:])
            exit(1)
        else:
            print("No such contributed binary:", argv[0], file=sys.stderr)
            exit(1)
    else:
        print("Usage: python3 -m pySecDecContrib --dirname | program argument ...", file=sys.stderr)
        exit(1)
