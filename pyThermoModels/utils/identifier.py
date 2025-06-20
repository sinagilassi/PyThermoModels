# import libs
from functools import wraps
# local


def add_attributes(**attrs):
    """Universal decorator for attaching attributes to any function."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        for key, value in attrs.items():
            setattr(wrapper, key, value)
        return wrapper
    return decorator
