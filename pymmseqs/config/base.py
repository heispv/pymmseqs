# pymmseqs/config/base.py

from abc import ABC, abstractmethod
from typing import Any, Dict, Union
from pathlib import Path
import yaml

from ..utils import (
    resolve_path,
    get_caller_dir,
    add_arg,
    add_twin_arg
)

class BaseConfig(ABC):

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        self._path_params = []
        self._required_files = []

    @abstractmethod
    def validate(self):
        """Validate the configuration parameters."""
        pass

    def to_dict(self, exclude_private: bool = True) -> Dict[str, Any]:
        """Convert config to dictionary, excluding None values.
        
        Args:
            exclude_private: If True, excludes attributes starting with '_' (like _defaults)
        """
        base_dict = {k: v for k, v in self.__dict__.items() 
                    if v is not None and (not exclude_private or not k.startswith('_'))}
        return base_dict

    @classmethod
    def from_yaml(cls, yaml_path: Union[str, Path]) -> 'BaseConfig':
        """
        Create a config instance from a YAML file.
        
        Args:
            yaml_path: Path to the YAML configuration file
            
        Returns:
            Config instance
            
        Raises:
            ValueError: If required fields are missing
            FileNotFoundError: If any input file doesn't exist
        """
        caller_dir = Path(get_caller_dir())
        yaml_path = resolve_path(yaml_path, caller_dir)
        
        with open(yaml_path, 'r') as f:
            config_dict = yaml.safe_load(f)
            
        return cls(**config_dict)
    
    def resolve_all_path(self, base_dir: Path) -> None:
        """
        Resolve all path specified in _path_params.
        
        Args:
            base_dir: Base directory for resolving relative path
        """
        for param in self._path_params:
            value = getattr(self, param, None)
            if value is None:
                continue

            if isinstance(value, list):
                resolved_values = [
                    str(resolve_path(path, base_dir))
                    for path in value
                ]
                setattr(self, param, resolved_values)
            else:
                resolved = str(resolve_path(value, base_dir))
                setattr(self, param, resolved)

    def validate_required_files(self) -> None:
        """
        Validate that all required files exist.
        
        Raises:
            FileNotFoundError: If any required file doesn't exist
        """
        for param in self._required_files:
            value = getattr(self, param, None)
            if value is None:
                continue

            if isinstance(value, list):
                for path in value:
                    if not Path(path).exists():
                        raise FileNotFoundError(f"Required file not found: {path}")
            else:
                if not Path(value).exists():
                    raise FileNotFoundError(f"Required file not found: {value}")

    def get_command_args(self, command_name: str) -> list:
        """
        Create command arguments based on the configuration.
        
        Returns:
            list: Command arguments starting with command name followed by parameters
        """
        # Create the command arguments starting with the command name from YAML
        args = [command_name]
        
        # Loop through all parameters and add the arguments
        for param_name, param_info in self._defaults.items():
            if param_info['required']:
                value = getattr(self, param_name)
                if isinstance(value, list):
                    args.extend(str(v) for v in value)
                else:
                    args.append(str(value))
            else:
                cmd_param = f"--{param_name.replace('_', '-')}"
                
                current_value = getattr(self, param_name)
                default_value = param_info['default']
                
                if param_info['twin']:
                    add_twin_arg(args, cmd_param, current_value, default_value, ',')
                else:
                    add_arg(args, cmd_param, current_value, default_value)
        
        return args