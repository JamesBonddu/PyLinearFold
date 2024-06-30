import re
import sys
from typing import Any, Dict, List, Optional, Type

class Value:
    def parse(self, text: str) -> None:
        raise NotImplementedError

    def parse_default(self) -> None:
        raise NotImplementedError

    def has_arg(self) -> bool:
        raise NotImplementedError

    def has_default(self) -> bool:
        raise NotImplementedError

    def is_container(self) -> bool:
        raise NotImplementedError

    def has_implicit(self) -> bool:
        raise NotImplementedError

    def get_default_value(self) -> str:
        raise NotImplementedError

    def get_implicit_value(self) -> str:
        raise NotImplementedError

    def default_value(self, value: str) -> 'Value':
        raise NotImplementedError

    def implicit_value(self, value: str) -> 'Value':
        raise NotImplementedError

class OptionException(Exception):
    def __init__(self, message: str):
        self.message = message

class OptionSpecException(OptionException):
    pass

class OptionParseException(OptionException):
    pass

class OptionExistsError(OptionSpecException):
    def __init__(self, option: str):
        super().__init__(f"Option ‘{option}’ already exists")

class InvalidOptionFormatError(OptionSpecException):
    def __init__(self, format: str):
        super().__init__(f"Invalid option format ‘{format}’")

class OptionNotExistsException(OptionParseException):
    def __init__(self, option: str):
        super().__init__(f"Option ‘{option}’ does not exist")

class MissingArgumentException(OptionParseException):
    def __init__(self, option: str):
        super().__init__(f"Option ‘{option}’ is missing an argument")

class OptionRequiresArgumentException(OptionParseException):
    def __init__(self, option: str):
        super().__init__(f"Option ‘{option}’ requires an argument")

class OptionNotHasArgumentException(OptionParseException):
    def __init__(self, option: str, arg: str):
        super().__init__(f"Option ‘{option}’ does not take an argument, but argument ‘{arg}’ given")

class OptionNotPresentException(OptionParseException):
    def __init__(self, option: str):
        super().__init__(f"Option ‘{option}’ not present")

class ArgumentIncorrectType(OptionParseException):
    def __init__(self, arg: str):
        super().__init__(f"Argument ‘{arg}’ failed to parse")

class StandardValue(Value):
    def __init__(self, value_type: Type, default: Optional[Any] = None):
        self.value_type = value_type
        self.default = default
        self.implicit = None
        self.store = None

    def parse(self, text: str) -> None:
        if text:
            self.store = self.value_type(text)
        elif self.implicit is not None:
            self.store = self.value_type(self.implicit)
        else:
            raise MissingArgumentException("Option requires an argument")

    def parse_default(self) -> None:
        if self.default is not None:
            self.store = self.value_type(self.default)

    def has_arg(self) -> bool:
        return True

    def has_default(self) -> bool:
        return self.default is not None

    def is_container(self) -> bool:
        return False

    def has_implicit(self) -> bool:
        return self.implicit is not None

    def get_default_value(self) -> str:
        return str(self.default)

    def get_implicit_value(self) -> str:
        return str(self.implicit)

    def default_value(self, value: str) -> 'Value':
        self.default = value
        return self

    def implicit_value(self, value: str) -> 'Value':
        self.implicit = value
        return self

class Options:
    def __init__(self, program: str, help_string: str = ""):
        self.program = program
        self.help_string = help_string
        self.options: Dict[str, StandardValue] = {}
        self.positional: List[str] = []
        self.next_positional = 0

    def add_options(self, group: str = "") -> 'OptionAdder':
        return OptionAdder(self, group)

    def add_option(self, group: str, s: str, l: str, desc: str, value: StandardValue, arg_help: str) -> None:
        if s:
            self.options[s] = value
        if l:
            self.options[l] = value

    def parse(self, argv: List[str]) -> None:
        i = 1
        while i < len(argv):
            arg = argv[i]
            if arg == "--":
                i += 1
                break

            match = re.match(r"--([a-zA-Z0-9][-_a-zA-Z0-9]*)(=(.*))?|-([a-zA-Z]+)", arg)
            if not match:
                if not self.consume_positional(arg):
                    argv[i] = arg
                    i += 1
                continue

            if match.group(4):
                short_options = match.group(4)
                for opt in short_options:
                    if opt in self.options:
                        if self.options[opt].has_arg():
                            if i + 1 < len(argv):
                                self.options[opt].parse(argv[i + 1])
                                i += 1
                            else:
                                raise MissingArgumentException(opt)
                        else:
                            self.options[opt].parse("")
                    else:
                        raise OptionNotExistsException(opt)
            elif match.group(1):
                long_option = match.group(1)
                if long_option in self.options:
                    if match.group(3):
                        self.options[long_option].parse(match.group(3))
                    elif self.options[long_option].has_arg():
                        if i + 1 < len(argv):
                            self.options[long_option].parse(argv[i + 1])
                            i += 1
                        else:
                            raise MissingArgumentException(long_option)
                    else:
                        self.options[long_option].parse("")
                else:
                    raise OptionNotExistsException(long_option)

            i += 1

        for opt in self.options.values():
            if opt.store is None and opt.has_default():
                opt.parse_default()

    def consume_positional(self, arg: str) -> bool:
        if self.next_positional < len(self.positional):
            opt = self.positional[self.next_positional]
            if opt in self.options:
                self.options[opt].parse(arg)
                self.next_positional += 1
                return True
        return False

    def parse_positional(self, options: List[str]) -> None:
        self.positional = options
        self.next_positional = 0

    def help(self) -> str:
        help_text = f"{self.help_string}\nUsage:\n  {self.program} [OPTION...]"
        if self.positional:
            help_text += " positional parameters"
        help_text += "\n\n"

        for opt in self.options.values():
            help_text += f"  --{opt.get_default_value()}\n"

        return help_text

class OptionAdder:
    def __init__(self, options: Options, group: str):
        self.options = options
        self.group = group

    def __call__(self, opts: str, desc: str, value: Optional[StandardValue] = None, arg_help: str = "") -> 'OptionAdder':
        match = re.match(r"(([a-zA-Z]),)?([a-zA-Z][-_a-zA-Z0-9]+)", opts)
        if not match:
            raise InvalidOptionFormatError(opts)

        s = match.group(2)
        l = match.group(3)

        if value is None:
            value = StandardValue(str)

        self.options.add_option(self.group, s, l, desc, value, arg_help)
        return self


if __name__ == "__main__":
    options = Options("program", "Help string")
    options.add_options()("o,option", "Option description", StandardValue(str).default_value("default"))
    options.parse_positional(["option"])
    options.parse(sys.argv)
    print(options.help())