from typing import Any, Dict, Optional, Type, TypeVar, cast

# Define a generic type variable
T = TypeVar('T')

class SingletonMeta(type):
    """
    This is a metaclass for creating a Singleton. By using this metaclass,
    we ensure that only one instance of the singleton class will be created.
    """
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            instance = super().__call__(*args, **kwargs)
            cls._instances[cls] = instance
        return cls._instances[cls]

class SingletonFactory(metaclass=SingletonMeta):
    """
    SingletonFactory itself is a singleton, ensuring only one instance exists.
    It can manage creation and access to other services or objects.
    """
    def __init__(self):
        if hasattr(self, '_initialized') and self._initialized:
            return
        self._initialized = True
        self._services = {} 

    def add_service(self, name, service):
        if name not in self._services:
            self._services[name] = service

    def get_service(self, name):
        return self._services.get(name)

    def remove_service(self, name):
        if name in self._services:
            del self._services[name]
            
    def get_service_as(self, name: str, cls: Type[T]) -> Optional[T]:
        service = self.get_service(name)
        if isinstance(service, cls):
            return cast(T, service)
        return None