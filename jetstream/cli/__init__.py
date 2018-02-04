""" This cli module is setup to allow some common arguments to always be
 accepted by the application, then dispatch the remaining arguments to
 subcommands. To expand the function of the command line utility, add another
 module to the "subcommands" package. It is required to have a function called
 "arg_parser" where the main parser will be passed so that you can add a
 subparser to it. Here is the bare minimum that the function needs:

```python
def arg_parser(subparser):
    parser = subparser.add_parser('SUBCOMMAND_NAME')
    parser.set_defaults(action=MAIN_FUNCTION)
```

Please, do not add a subcommand name that conflicts with an existing subcommand.

When the cli is invoked with the given SUBCOMMAND_NAME, the function given in
MAIN_FUNCTION will be passed a argparse namespace object, containing the main
arguments as well as subcommand arguments. Subcommand arguments are independent
from other subcommands.

 """