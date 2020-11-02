from functools import wraps, partial


class ApproximationError(RuntimeError):
    pass


class Model:
    def __init__(self):
        self.approximations = set()

    def with_approximation(self, name):
        self.approximations.add(name)
        return self


class _Wrapper:
    def __init__(self, func):
        self.func = func
        self.approximations = set()
        self.original_doc = func.__doc__

    def __get__(self, obj, cls=None):
        @wraps(self.func)
        def call(*args, **kwargs):
            if not self.approximations.issubset(obj.approximations):
                raise ApproximationError
            return self.func(obj, *args, **kwargs)

        return call

    def update_doc(self):
        self.func.__doc__ = f'**Required approximations:** {", ".join(self.approximations)}\n{self.original_doc}'


class approximation:  # noqa
    def __init__(self, name):
        self.approx_name = name

    def __call__(self, func):
        wrapper = func if isinstance(func, _Wrapper) else _Wrapper(func)
        wrapper.approximations.add(self.approx_name)
        wrapper.update_doc()
        return wrapper
