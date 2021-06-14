if __name__ == "__main__":
    import sys
    from . import dirname
    argv = sys.argv[1:]
    if argv == ["--dirname"]:
        print(dirname)
    else:
        print("Usage: python3 -m pySecDecContrib --dirname", file=sys.stderr)
        exit(1)
