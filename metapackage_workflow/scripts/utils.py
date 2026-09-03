import os
import hashlib
import pickle
import gzip
import shutil
import joblib
import logging

from contextlib import contextmanager

from glob import glob
from typing import Callable, List, Union

from functools import wraps


# Source: https://stackoverflow.com/questions/24983493/tracking-progress-of-joblib-parallel-execution
@contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager which provides progress bar for joblib operations."""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


class CacheManager:
    def __init__(self, logger: logging.Logger = None, cache_dir: str = 'cache'):
        self.cache_dir = cache_dir
        self._script_cached_files = []
        
        # Use provided logger or create a default one
        if logger is None:
            self.logger = logging.getLogger(__name__)
            # Only configure if it doesn't have handlers (avoid duplicate handlers)
            if not self.logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
                self.logger.setLevel(logging.INFO)
        else:
            self.logger = logger

        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
            self.logger.debug(f"Created cache directory: {self.cache_dir}")

    def cache_results(self, func):
        def wrapper(*args, **kwargs):
            cache_key = f"{func.__name__}_{hashlib.md5(str(args).encode() + str(kwargs).encode()).hexdigest()}.pkl"
            cache_path = os.path.join(self.cache_dir, cache_key)

            if cache_path not in self._script_cached_files:
                self._script_cached_files.append(cache_path)

            if os.path.exists(cache_path):
                with open(cache_path, 'rb') as cache_file:
                    self.logger.debug(f"Loading cached result from {cache_path}")
                    return pickle.load(cache_file)

            result = func(*args, **kwargs)
            with open(cache_path, 'wb') as cache_file:
                pickle.dump(result, cache_file)
                self.logger.debug(f"Saving result to cache at {cache_path}")

            return result
        return wrapper
    

    def _cleanup_cached_files(self):
        self.logger.debug(f"[{self.cache_dir}] Attempting to clean up cached results for this script...")
        for path in self._script_cached_files:
            try:
                if os.path.exists(path):
                    os.remove(path)
                    self.logger.debug(f"[{self.cache_dir}] Deleted cached file: {path}")
                else:
                    self.logger.debug(f"[{self.cache_dir}] Cached file not found (already deleted or never created): {path}")
            except Exception as e:
                self.logger.error(f"[{self.cache_dir}] Error deleting {path}: {e}")
        self._script_cached_files.clear()



def handle_gzip(logger: logging.Logger) -> Callable:
    """Decorator to handle gzip files in all function arguments (both positional and keyword)."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            temp_files = []  # Track all temporary files we create
            modified_args = list(args)  # Convert to list so we can modify
            
            try:
                # Check positional arguments for .gz files
                for i, arg in enumerate(args):
                    if isinstance(arg, str) and arg.endswith('.gz'):
                        temp_path = _create_temp_file(arg, logger)
                        temp_files.append(temp_path)
                        modified_args[i] = temp_path
                        logger.debug(f"Updated positional arg {i} to use temporary file: {temp_path}")
                
                # Check keyword arguments for .gz files
                for arg_name, file_path in kwargs.items():
                    if isinstance(file_path, str) and file_path.endswith('.gz'):
                        temp_path = _create_temp_file(file_path, logger)
                        temp_files.append(temp_path)
                        kwargs[arg_name] = temp_path
                        logger.debug(f"Updated keyword arg '{arg_name}' to use temporary file: {temp_path}")

                # Execute the original function with modified arguments
                result = func(*modified_args, **kwargs)
                
            except Exception as e:
                logger.error(f"Error in function {func.__name__}: {e}")
                raise
            finally:
                # Clean up all temporary files
                _cleanup_temp_files(temp_files, logger)

            return result
        return wrapper
    return decorator

def _create_temp_file(gz_file_path: str, logger: logging.Logger) -> str:
    """Helper function to create and decompress a temporary file."""
    import tempfile
    
    # Get the original extension
    original_name = os.path.basename(gz_file_path)
    if original_name.endswith('.gz'):
        uncompressed_name = original_name[:-3]  # Remove .gz
    else:
        uncompressed_name = original_name
    
    # Create unique temporary file with original extension
    with tempfile.NamedTemporaryFile(
        suffix=f"_{uncompressed_name}",
        delete=False,
        prefix="temp_gzip_"
    ) as temp_file:
        temp_path = temp_file.name
    
    # Decompress to temporary file
    logger.debug(f"Decompressing {gz_file_path} to {temp_path}")
    try:
        with gzip.open(gz_file_path, 'rb') as f_in:
            with open(temp_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    except Exception as e:
        logger.error(f"Error decompressing '{gz_file_path}': {e}")
        # Clean up failed temp file
        if os.path.exists(temp_path):
            os.remove(temp_path)
        raise
    
    return temp_path

def _cleanup_temp_files(temp_files: list, logger: logging.Logger):
    """Helper function to clean up temporary files."""
    for temp_path in temp_files:
        if os.path.exists(temp_path):
            try:
                os.remove(temp_path)
                logger.debug(f"Cleaned up temporary file: {temp_path}")
            except OSError as e:
                logger.warning(f"Could not remove temporary file {temp_path}: {e}")

def glob_genomes(
    genomes_dir: str,
) -> Union[List[str], str]:

    genomes_files = glob(genomes_dir, recursive=True, include_hidden=True)
    genomes_files = sorted(genomes_files, key=lambda x: x.split('/')[-1])
    
    genomes_files = [f for f in genomes_files if os.path.isfile(f)]

    if len(genomes_files) == 1:
        genomes_files = genomes_files[0]  # If only one file, return it directly
    
    return genomes_files



@contextmanager
def named_temp_directory(name, dir=None):
    """
    Context manager that creates a temporary directory with a predictable name.
    
    Parameters:
    -----------
    prefix : str
        Prefix for the directory name
    suffix : str  
        Suffix for the directory name
    dir : str
        Parent directory (defaults to system temp)
    
    Yields:
    -------
    str
        Path to the temporary directory
    """
    if dir:
        temp_dir = os.path.join(dir, name)
    else:
        temp_dir = name

    os.makedirs(temp_dir, exist_ok=True)

    try:
        yield temp_dir
    finally:
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)