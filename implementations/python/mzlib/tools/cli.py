import sys
import traceback
import click

from mzlib.spectrum_library import SpectrumLibrary
from mzlib.index import MemoryIndex, SQLIndex
from mzlib.backends.text import TextSpectralLibraryWriter

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    '''A collection of utilities for inspecting and manipulating
    spectral libraries.
    '''
    pass


@main.command("describe", short_help=("Produce a minimal textual description"
                                     " of a spectral library"))
@click.argument('path', type=click.Path(exists=True))
@click.option("-d", "--diagnostics", is_flag=True,
              help="Run more diagnostics, greatly increasing runtime but producing additional information")
def describe(path, diagnostics=False):
    '''Produces a minimal textual description of a spectral library.
    '''
    click.echo("Describing \"%s\"" % (path,))
    if SQLIndex.exists(path):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex
    library = SpectrumLibrary(filename=path, index_type=index_type)
    click.echo(f"Format: {library.format}")
    click.echo(f"Size: {library.__len__()}")
    fh = click.open_file("-", 'wt')
    TextSpectralLibraryWriter(fh).write_header(library.backend)



@main.command("convert", short_help=("Convert a spectral library from one format to another"))
@click.argument('inpath', type=click.Path(exists=True))
@click.argument("outpath", type=click.Path())
@click.option("-f", "--format", type=click.Choice(["text", "json"]), default="text")
def convert(inpath, outpath, format=None):
    '''Convert a spectral library from one format to another. If `outpath` is `-`,
    instead of writing to file, data will instead be sent to STDOUT.
    '''
    if format is None:
        format = "text"
    if SQLIndex.exists(inpath):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex
    library = SpectrumLibrary(filename=inpath, index_type=index_type)
    fh = click.open_file(outpath, mode='w')
    library.write(fh, format)


if __name__ == "__main__":
    main()
