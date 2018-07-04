#!/usr/bin/env python
# encoding: utf-8
#
# logger.py
#
# Created by José Sánchez-Gallego on 17 Sep 2017.


from __future__ import absolute_import, division, print_function

import datetime
import logging
import os
import pathlib
import re
import shutil
import sys
import traceback
import warnings
from logging import PercentStyle
from logging.handlers import TimedRotatingFileHandler

import click
from pygments import highlight
from pygments.formatters import TerminalFormatter
from pygments.lexers import get_lexer_by_name

from ..exceptions import GohanWarning


# from textwrap import TextWrapper


# Adds custom log level for print and twisted messages
PRINT = 15
logging.addLevelName(PRINT, 'PRINT')


def print_log_level(self, message, *args, **kws):
    self._log(PRINT, message, args, **kws)


logging.Logger._print = print_log_level


# Adds custom log level for important messages
IMPORTANT = 25
logging.addLevelName(IMPORTANT, 'IMPORTANT')


def important(self, message, *args, **kws):
    self._log(IMPORTANT, message, args, **kws)


logging.Logger.important = important


def print_exception_formatted(type, value, tb):
    """A custom hook for printing tracebacks with colours."""

    tbtext = ''.join(traceback.format_exception(type, value, tb))
    lexer = get_lexer_by_name('pytb', stripall=True)
    formatter = TerminalFormatter()
    sys.stderr.write(highlight(tbtext, lexer, formatter))


def colored_formatter(record):
    """Prints log messages with colours."""

    colours = {'info': ('blue', 'normal'),
               'debug': ('magenta', 'normal'),
               'warning': ('yellow', 'normal'),
               'print': ('green', 'normal'),
               'error': ('red', 'bold'),
               'important': ('green', 'bold')}

    levelname = record.levelname.lower()

    if levelname == 'error':
        return

    if levelname.lower() in colours:
        levelname_color = colours[levelname][0]
        bold = True if colours[levelname][1] == 'bold' else False
        header = click.style('[{}]: '.format(levelname.upper()), levelname_color, bold=bold)

    message = record.getMessage()

    warning_category = re.match(r'^(\w+Warning:)(.*)', message)
    if warning_category is not None:
        warning_category_colour = click.style(warning_category.groups()[0], 'cyan')
        message = warning_category_colour + warning_category.groups()[1]

    # sub_level = re.match('(\[.+\]:)(.*)', message)
    # if sub_level is not None:
    #     sub_level_name = click.style(sub_level.groups()[0], 'red')
    #     message = '{}{}'.format(sub_level_name, ''.join(sub_level.groups()[1:]))

    # if len(message) > 79:
    #     tw = TextWrapper()
    #     tw.width = 79
    #     tw.subsequent_indent = ' ' * (len(record.levelname) + 2)
    #     tw.break_on_hyphens = False
    #     message = '\n'.join(tw.wrap(message))

    sys.__stdout__.write('{}{}\n'.format(header, message))
    sys.__stdout__.flush()

    return


class MyFormatter(logging.Formatter):

    warning_fmt = '%(asctime)s - %(levelname)s: %(message)s [%(origin)s]'
    info_fmt = '%(asctime)s - %(levelname)s - %(message)s [%(funcName)s @ ' + \
        '%(filename)s]'

    ansi_escape = re.compile(r'\x1b[^m]*m')

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        # format_orig = self._fmt

        # Replace the original format with one customized by logging level

        if record.levelno == logging.DEBUG:
            self._style = PercentStyle(MyFormatter.info_fmt)

        elif record.levelno == logging.getLevelName('PRINT'):
            self._style = PercentStyle(MyFormatter.info_fmt)

        elif record.levelno == logging.getLevelName('IMPORTANT'):
            self._style = PercentStyle(MyFormatter.info_fmt)

        elif record.levelno == logging.INFO:
            self._style = PercentStyle(MyFormatter.info_fmt)

        elif record.levelno == logging.ERROR:
            self._style = PercentStyle(MyFormatter.info_fmt)

        elif record.levelno == logging.WARNING:
            self._style = PercentStyle(MyFormatter.warning_fmt)

        record.msg = self.ansi_escape.sub('', record.msg)

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        # self._fmt = format_orig

        return result


Logger = logging.getLoggerClass()
fmt = MyFormatter()


class LoggerStdout(object):
    """A pipe for stdout to a logger."""

    def __init__(self, level):
        self.level = level

    def write(self, message):

        if message != '\n':
            self.level(message)

    def flush(self):
        pass


class MyLogger(Logger):
    """This class is used to set up the logging system.

    The main functionality added by this class over the built-in
    logging.Logger class is the ability to keep track of the origin of the
    messages, the ability to enable logging of warnings.warn calls and
    exceptions, and the addition of colorized output and context managers to
    easily capture messages to a file or list.

    """

    INFO = 15
    IMPORTANT = 25

    # The default actor to log to. It is set by the set_actor() method.
    _actor = None

    def save_log(self, path):
        shutil.copyfile(self.log_filename, os.path.expanduser(path))

    def _warn(self, message, category=None, stacklevel=1):
        """Overrides `warnings.warn`

        Before calling the original `warnings.warn` function it makes sure
        the warning is redirected to the correct ``showwarning`` function.

        """

        if issubclass(category, GohanWarning):
            warnings.showwarning = self._show_warning
        else:
            warnings.showwarning = warnings._showwarning_orig

        warnings._original_warn(message, category=category, stacklevel=stacklevel)

    def _show_warning(self, *args, **kwargs):

        warning = args[0]

        message = '{0}: {1}'.format(warning.__class__.__name__, args[0])
        mod_path = args[2]

        mod_name = None
        mod_path, ext = os.path.splitext(mod_path)
        for name, mod in sys.modules.items():
            mod_file = getattr(mod, '__file__', '')
            if mod_file is not None:
                path = os.path.splitext(mod_file)[0]
                if path == mod_path:
                    mod_name = mod.__name__
                    break

        if mod_name is not None:
            warning.logger.warning(message, extra={'origin': mod_name})
        else:
            warning.logger.warning(message, extra={'origin': 'no_module'})

    def _catch_exceptions(self, exctype, value, tb):
        """Catches all exceptions and logs them."""

        # Now we log it.
        self.error('Uncaught exception', exc_info=(exctype, value, tb))

        # First, we print to stdout with some colouring.
        print_exception_formatted(exctype, value, tb)

    def _set_defaults(self, log_level=logging.INFO, redirect_stdout=False):
        """Reset logger to its initial state."""

        # Remove all previous handlers
        for handler in self.handlers[:]:
            self.removeHandler(handler)

        # Set levels
        self.setLevel(logging.DEBUG)

        # Set up the stdout handler
        self.fh = None
        self.sh = logging.StreamHandler()
        self.sh.emit = colored_formatter
        self.addHandler(self.sh)

        self.sh.setLevel(log_level)

        self.enable_warnings()

        # Redirects all stdout to the logger
        if redirect_stdout:
            sys.stdout = LoggerStdout(self._print)

        # Catches exceptions
        sys.excepthook = self._catch_exceptions

    def enable_warnings(self):
        """Redirects warnings to the log."""

        warnings.showwarning = self._show_warning

        warnings._original_warn = warnings.warn
        warnings.warn = self._warn

    def disable_warnings(self):
        """Restores normal warning system."""

        warnings.showwarning = warnings._show_warning
        warnings.warn = warnings._original_warn

    def start_file_logger(self, name, log_file_level=logging.DEBUG, log_file_path='~/'):
        """Start file logging."""

        log_file_path = pathlib.Path(log_file_path).expanduser() / '{}.log'.format(name)
        logdir = log_file_path.parent

        try:
            logdir.mkdir(parents=True, exist_ok=True)

            # If the log file exists, backs it up before creating a new file handler
            if log_file_path.exists():
                strtime = datetime.datetime.utcnow().strftime('%Y-%m-%d_%H:%M:%S')
                shutil.move(str(log_file_path), str(log_file_path) + '.' + strtime)

            self.fh = TimedRotatingFileHandler(str(log_file_path), when='midnight', utc=True)
            self.fh.suffix = '%Y-%m-%d_%H:%M:%S'
        except (IOError, OSError) as ee:
            warnings.warn('log file {0!r} could not be opened for writing: '
                          '{1}'.format(log_file_path, ee), RuntimeWarning)
        else:
            self.fh.setFormatter(fmt)
            self.addHandler(self.fh)
            self.fh.setLevel(log_file_level)

        self.log_filename = log_file_path

    def saveLog(self, path):
        shutil.copyfile(self.log_filename, os.path.expanduser(path))


logging.setLoggerClass(MyLogger)
log = logging.getLogger(__name__)
log._set_defaults()  # Inits sh handler
