"""The cli_front is a the command line utility that is used to list all the
   accessable command line utilities and to call the command line utility
   you want to run."""
import sys, argparse, pkgutil
from importlib import import_module
import os.path
import seqtools.cli.utilities

def main():
  # Get the tasks available
  cache_argv = sys.argv
  front_end_args = sys.argv[:2]
  back_end_args = sys.argv[1:]
  # only use the front end arguments
  sys.argv = front_end_args

  args = do_args()
  sys.argv = cache_argv #put back the arguments

  # Utilities in the utility directory that will be available must have
  #    external_cmd(string) available
  task_module = import_module('seqtools.cli.utilities.'+args.task)
  task_module.external_cmd(" ".join(back_end_args))

def do_args ():
  util_path = os.path.dirname(seqtools.cli.utilities.__file__)
  packlist = [name for _, name, _ in pkgutil.iter_modules([util_path])]
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('task',choices=packlist,help="Specify which task to execute")
  args = parser.parse_args()
  return args

if __name__=="__main__":
  main()
