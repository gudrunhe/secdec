version = __version__ = '1.2'
authors = __authors__ = 'Sophia Borowka, Gudrun Heinrich, Stephan Jahn, Stephen Jones, Matthias Kerner, Johannes Schlenk, Tom Zirke'

from subprocess import check_output as shell_call
import os
# The following line is replaced by `setup.py` --> hard-code the commit id in distributions
git_id = shell_call(['git','rev-parse','HEAD'], cwd=os.path.split(os.path.abspath(__file__))[0]).strip()
# in python3, the return type of `shell_call` may be `bytes` but we need `str`
if not isinstance(git_id, str):
    git_id = git_id.decode()
