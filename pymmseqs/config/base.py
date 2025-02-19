# pymmseqs/config/base.py

from abc import ABC, abstractmethod
from typing import Any, Dict, Union
from pathlib import Path
import yaml

from ..utils import (
    resolve_path,
    get_caller_dir,
    add_arg
)

class BaseConfig(ABC):

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        self._defaults = {}

    @abstractmethod
    def _validate(self):
        """
        Validate the configuration parameters.
        
        This method must be implemented by subclasses to define their specific
        validation rules.
        
        Raises:
            ValueError: If any parameter fails validation
        """
        pass

    def to_dict(self, exclude_private: bool = True) -> Dict[str, Any]:
        """
        Convert configuration to a dictionary, excluding None values.
        
        Args:
            exclude_private (bool): If True, excludes attributes starting with '_'
                                  like _defaults and other internal attributes
        
        Returns:
            Dict[str, Any]: Dictionary containing the configuration parameters
        """
        base_dict = {k: v for k, v in self.__dict__.items() 
                    if v is not None and (not exclude_private or not k.startswith('_'))}
        return base_dict

    @classmethod
    def from_yaml(cls, yaml_path: Union[str, Path]) -> 'BaseConfig':
        """
        Create a configuration instance from a YAML file.
        
        Args:
            yaml_path (Union[str, Path]): Path to the YAML configuration file
            
        Returns:
            BaseConfig: New instance of the configuration class
            
        Raises:
            ValueError: If required fields are missing
            FileNotFoundError: If the YAML file or any required input files don't exist
        """
        caller_dir = Path(get_caller_dir())
        yaml_path = resolve_path(yaml_path, caller_dir)
        
        with open(yaml_path, 'r') as f:
            config_dict = yaml.safe_load(f)
            
        return cls(**config_dict)
    
    def _resolve_all_path(self, base_dir: Path) -> None:
        """
        Resolve all paths specified in _defaults relative to a base directory.
        
        This method handles both single paths and lists of paths, converting them
        to absolute paths based on the provided base directory.
        
        Args:
            base_dir (Path): Base directory for resolving relative paths
        """
        # Resolve all paths using _defaults
        for param_name, param_info in self._defaults.items():
            if param_info['type'] == 'path':
                value = getattr(self, param_name)
                if value:
                    if isinstance(value, list):
                        resolved_values = [str(resolve_path(v, base_dir)) for v in value]
                        setattr(self, param_name, resolved_values)
                    else:
                        resolved = str(resolve_path(value, base_dir))
                        setattr(self, param_name, resolved)

    def _check_required_files(self) -> None:
        """
        Verify that all required files exist.
        
        Checks both single files and lists of files marked as required in _defaults.
        
        Raises:
            ValueError: If a required file parameter is not set
            FileNotFoundError: If any required file doesn't exist
        """
        for param_name, param_info in self._defaults.items():
            if param_info['required'] and param_info['should_exist']:
                value = getattr(self, param_name)
                if value is None:
                    raise ValueError(f"Required file is not set: {param_name}")

            if isinstance(value, list):
                for path in value:
                    if not Path(path).exists():
                        raise FileNotFoundError(f"Required file not found: {path}")
            else:
                if not Path(value).exists():
                    raise FileNotFoundError(f"Required file not found: {value}")

    def _validate_choices(self):
        """
        Validate parameters against their allowed choices.
        
        Checks all parameters that have defined choices in _defaults to ensure
        they contain valid values. Skips optional parameters that are set to
        their default values.
        
        Raises:
            ValueError: If any parameter has an invalid value
        """
        for param_name, param_info in self._defaults.items():
            value = getattr(self, param_name)
            
            # Skip optional parameters with default values
            if not param_info['required'] and value == param_info['default']:
                continue
            
            # Validate choices if they exist
            if param_info['choices'] is not None:
                if value not in param_info['choices']:
                    raise ValueError(
                        f"{param_name} is {value} but must be one of {param_info['choices']}"
                    )

    def _get_command_args(self, command_name: str) -> list:
        """
        Generate command-line arguments for MMseqs2 execution.
        
        Creates a list of command-line arguments based on the configuration,
        handling different parameter types (boolean, twin, comma-separated strings)
        appropriately.

        - Check the related .yaml file for more information about the _defaults parameter for each command.
        
        Args:
            command_name (str): Name of the MMseqs2 command to execute
            
        Returns:
            list: Command arguments starting with command name followed by parameters
        """
        # Create the command arguments starting with the command name from YAML
        args = [command_name.replace('_', '-')]
        
        # Loop through all parameters and add the arguments
        for param_name, param_info in self._defaults.items():
            if param_info['required']:
                value = getattr(self, param_name)
                if isinstance(value, list):
                    args.append(','.join(map(str, value)))
                else:
                    args.append(str(value))
            else:
                # Get current and default values
                current_value = getattr(self, param_name)
                default_value = param_info['default']
                
                # Create parameter flag (handle single character parameters differently)
                cmd_param = f"-{param_name}" if len(param_name) == 1 else f"--{param_name.replace('_', '-')}"
                
                # Handle different parameter types
                if param_info['twin']:
                    # For twin parameters, compare as strings
                    if str(current_value) != str(default_value):
                        add_arg(args, cmd_param, current_value, default_value)
                        
                elif param_info['type'] == "comma_separated_str":
                    # For comma-separated strings, compare as lists
                    current_list = [item.strip() for item in str(current_value).split(",")]
                    default_list = [item.strip() for item in str(default_value).split(",")]
                    if current_list != default_list:
                        cleaned_value = ",".join(current_list)
                        args.extend([cmd_param, cleaned_value])
                            
                elif isinstance(current_value, bool):
                    # For boolean values, only add if different from default
                    if current_value != default_value:
                        args.extend([cmd_param, "1" if current_value else "0"])
                        
                else:
                    # For all other types, compare directly
                    if current_value != default_value:
                        args.extend([cmd_param, str(current_value)])
        
        return args
