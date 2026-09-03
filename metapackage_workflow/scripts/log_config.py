# log_config.py

import os
import logging

def configure_logger(file_path=None):
    if file_path is not None:
        logger_name = os.path.splitext(os.path.basename(file_path))[0]
    else:
        logger_name = __name__
    logger = logging.getLogger(logger_name)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s [%(name)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    return logger