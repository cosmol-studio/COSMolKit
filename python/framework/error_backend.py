"""The create_exception backend; no chemistry or business error inventory."""
import json
import types


def install(exception, payload):
    def no_public_constructor(cls, *args, **kwargs):
        raise TypeError(f'cannot create {cls.__name__} instances directly')

    def getter(descriptor):
        def read(self):
            return descriptor.__get__(self._ck_error_payload, payload)
        read.__wrapped__ = descriptor
        return property(read, doc=descriptor.__doc__)

    def method(descriptor):
        def call(self, *args, **kwargs):
            return descriptor.__get__(self._ck_error_payload, payload)(*args, **kwargs)
        call.__doc__ = descriptor.__doc__
        return call

    def source(function):
        def read(self):
            return function(self)
        read.__wrapped__ = function
        return property(read, doc=function.__doc__)

    for name, descriptor in vars(payload).items():
        if name.startswith('_') and name != '__str__':
            continue
        if isinstance(descriptor, types.GetSetDescriptorType):
            setattr(exception, name, getter(descriptor))
        elif isinstance(descriptor, (types.MethodDescriptorType, types.WrapperDescriptorType)):
            setattr(exception, name, method(descriptor))
        elif isinstance(descriptor, staticmethod):
            function = descriptor.__func__
            doc = (function.__doc__ or '').strip()
            if doc.startswith('{') and json.loads(doc).get('role') in {'error_source', 'error_payload'}:
                setattr(exception, name, source(function))
            else:
                raise TypeError(f'Unsupported error projection static method: {name}')
    exception.__name__ = payload.__name__
    exception.__qualname__ = payload.__qualname__
    exception.__module__ = payload.__module__
    exception.__doc__ = payload.__doc__
    exception.__new__ = staticmethod(no_public_constructor)
