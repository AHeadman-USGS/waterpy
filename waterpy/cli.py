"""Main `waterpy` command line interface"""


import click
import sys, warnings
import traceback

from waterpy.main import waterpy


class Options:
    def __init__(self):
        self.verbose = False
        self.show = False


# Create a decorator to pass options to each command
pass_options = click.make_pass_decorator(Options, ensure=True)


@click.group()
@click.option("-v", "--verbose", is_flag=True,
              help="Print model run details.")
@click.option("-s", "--show", is_flag=True,
              help="Show output plots.")
@click.pass_context
def main(options, verbose, show):
    """waterpy is a command line tool for a rainfall-runoff
    model that predicts the amount of water flow in rivers.
    """
    options.verbose = verbose
    options.show = show
    warnings.filterwarnings("ignore")


@main.command()
@click.argument("configfile", type=click.Path(exists=True))
@pass_options
def run(options, configfile):
    """Run waterpy with a model configuration file.

    The model configuration file contains the specifications for a model run.
    This command takes in the path to model configuration file.
    """
    try:
        click.echo("Running model...")
        waterpy(configfile, options)
        click.echo("Finished!")
        click.echo("Output saved as specified in the model config file.")
    except Exception as err:
        click.echo(err, traceback.print_exc())
        sys.exit(1)

    if options.verbose:
        click.echo("Verbose on")
    if options.show:
        click.echo("Show on")
