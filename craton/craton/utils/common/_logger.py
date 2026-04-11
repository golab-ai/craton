import copy
import logging
import logging.handlers
from pathlib import Path

# LOG_DIR = Path(__file__).parent.parent / "logs"
# Path(LOG_DIR).mkdir(parents=True, exist_ok=True)
LOG_DIR = Path("/tmp")


def _get_logger(fname: str, log_directory: str = None) -> logging.Logger:
    logger = logging.getLogger(fname)
    logger.propagate = False
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)s %(levelname)s %(filename)s:%(lineno)d | %(message)s")
    colored_formatter = ColoredFormatter("%(asctime)s %(levelname)s %(filename)s:%(lineno)d | %(message)s")

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(colored_formatter)
    logger.addHandler(stream_handler)

    if log_directory is None:
        logfile = str(LOG_DIR / fname)
    else:
        logfile = str(Path(log_directory) / fname)
    try:
        file_handler = logging.handlers.TimedRotatingFileHandler(filename=logfile, when="d", backupCount=90)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    except:
        logger.warning(f"Cannot open {logfile} for logging")

    return logger


BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

# The background is set with 40 plus the number of the color, and the foreground with 30

# These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

COLORS = {"WARNING": YELLOW, "INFO": WHITE, "DEBUG": BLUE, "CRITICAL": YELLOW, "ERROR": RED}


class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color=True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, rec):
        record = copy.copy(rec)
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
            return logging.Formatter.format(self, record)



def reset_logger(logger, dirname="current", fname=None, log_level="DEBUG"):
    # Reset root logger setting
    logger.setLevel(log_level)

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Set stream handler
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(filename)s:%(lineno)d | %(message)s")
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if fname is None:
        return

    # Set locations for log files
    if dirname == "current":
        log_dir = Path.cwd()
    else:
        logger.warning(f"Directory for logging files {dirname} not implemented, using current")
        log_dir = Path.cwd()

    # Set file handler
    logfile = str(log_dir / fname)
    try:
        file_handler = logging.handlers.TimedRotatingFileHandler(filename=logfile, when="d", backupCount=90)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    except:
        logger.warning(f"Cannot open {logfile} for logging")
