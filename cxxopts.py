import re
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

class OptionException(Exception):
    def __init__(self, message: str):
        self.message = message

    def __str__(self):
        return self.message

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

class StandardValue(Value):
    def __init__(self, store: Optional[Any] = None):
        self.store = store if store is not None else []
        self.default = False
        self.default_value_str = ""
        self.implicit = False
        self.implicit_value_str = ""

    def parse(self, text: str) -> None:
        if self.implicit and not text:
            self._parse_value(self.implicit_value_str)
        else:
            self._parse_value(text)

    def _parse_value(self, text: str) -> None:
        try:
            if isinstance(self.store, bool):
                self.store = True
            elif isinstance(self.store, list):
                self.store.append(self._convert(text))
            else:
                self.store = self._convert(text)
        except ValueError:
            raise ArgumentIncorrectType(text)

    def _convert(self, text: str) -> Any:
        if isinstance(self.store, bool):
            return True
        elif isinstance(self.store, list):
            return type(self.store[0])(text)
        else:
            return type(self.store)(text)

    def is_container(self) -> bool:
        return isinstance(self.store, list)

    def parse_default(self) -> None:
        self._parse_value(self.default_value_str)

    def has_arg(self) -> bool:
        return not isinstance(self.store, bool)

    def has_default(self) -> bool:
        return self.default

    def has_implicit(self) -> bool:
        return self.implicit

    def get_default_value(self) -> str:
        return self.default_value_str

    def get_implicit_value(self) -> str:
        return self.implicit_value_str

    def default_value(self, value: str) -> 'Value':
        self.default = True
        self.default_value_str = value
        return self

    def implicit_value(self, value: str) -> 'Value':
        self.implicit = True
        self.implicit_value_str = value
        return self

class OptionDetails:
    def __init__(self, description: str, value: Value):
        self.description = description
        self.value = value
        self.count = 0

    def parse(self, text: str) -> None:
        self.value.parse(text)
        self.count += 1

    def parse_default(self) -> None:
        self.value.parse_default()

    def has_arg(self) -> bool:
        return self.value.has_arg()

    def get_count(self) -> int:
        return self.count

class Options:
    def __init__(self, program: str, help_string: str = ""):
        self.program = program
        self.help_string = help_string
        self.options: Dict[str, OptionDetails] = {}
        self.positional: List[str] = []
        self.next_positional = 0
        self.positional_set: set[str] = set()
        self.help: Dict[str, Dict[str, Any]] = {}

    def parse(self, args: List[str]) -> None:
        current = 0
        while current < len(args):
            arg = args[current]
            if arg == "--":
                current += 1
                break

            match = re.match(r"--(\w[\w-]*)(=(.*))?|-(\w+)", arg)
            if not match:
                if not self._consume_positional(arg):
                    args[current] = arg
                    current += 1
                continue

            if match.group(4):
                short_options = match.group(4)
                for opt in short_options:
                    self._parse_short_option(opt, args, current)
            elif match.group(1):
                long_option = match.group(1)
                value = match.group(3)
                self._parse_long_option(long_option, value, args, current)

            current += 1

        for opt in self.options.values():
            if not opt.get_count() and opt.value.has_default():
                opt.parse_default()

    def _parse_short_option(self, opt: str, args: List[str], current: int) -> None:
        option_details = self.options.get(opt)
        if not option_details:
            raise OptionNotExistsException(opt)

        if option_details.has_arg():
            if current + 1 < len(args) and not args[current + 1].startswith("-"):
                option_details.parse(args[current + 1])
                current += 1
            elif option_details.value.has_implicit():
                option_details.parse("")
            else:
                raise MissingArgumentException(opt)
        else:
            option_details.parse("")

    def _parse_long_option(self, opt: str, value: Optional[str], args: List[str], current: int) -> None:
        option_details = self.options.get(opt)
        if not option_details:
            raise OptionNotExistsException(opt)

        if value is not None:
            if not option_details.has_arg():
                raise OptionNotHasArgumentException(opt, value)
            option_details.parse(value)
        else:
            if option_details.has_arg():
                if current + 1 < len(args) and not args[current + 1].startswith("-"):
                    option_details.parse(args[current + 1])
                    current += 1
                elif option_details.value.has_implicit():
                    option_details.parse("")
                else:
                    raise MissingArgumentException(opt)
            else:
                option_details.parse("")

    def _consume_positional(self, arg: str) -> bool:
        while self.next_positional < len(self.positional):
            opt = self.positional[self.next_positional]
            option_details = self.options.get(opt)
            if option_details:
                if not option_details.value.is_container() and option_details.get_count() == 0:
                    option_details.parse(arg)
                    self.next_positional += 1
                    return True
                elif option_details.value.is_container():
                    option_details.parse(arg)
                    return True
            self.next_positional += 1
        return False

    def add_options(self, group: str = "") -> 'OptionAdder':
        return OptionAdder(self, group)

    def add_option(self, group: str, short: str, long: str, desc: str, value: Value, arg_help: str) -> None:
        option_details = OptionDetails(desc, value)
        if short:
            self.options[short] = option_details
        if long:
            self.options[long] = option_details

        if group not in self.help:
            self.help[group] = {"options": []}
        self.help[group]["options"].append({
            "short": short,
            "long": long,
            "desc": desc,
            "has_arg": value.has_arg(),
            "has_default": value.has_default(),
            "default_value": value.get_default_value(),
            "has_implicit": value.has_implicit(),
            "implicit_value": value.get_implicit_value(),
            "arg_help": arg_help,
            "is_container": value.is_container()
        })

    def parse_positional(self, options: List[str]) -> None:
        self.positional = options
        self.next_positional = 0
        self.positional_set.update(options)

    def help_one_group(self, group: str) -> str:
        if group not in self.help:
            return ""

        options = self.help[group]["options"]
        longest = max(len(f"-{opt['short']}, --{opt['long']}" if opt['short'] else f"--{opt['long']}") for opt in options)
        longest = min(longest, 30)
        result = f" {group} options:\n" if group else ""

        for opt in options:
            if opt["is_container"] and opt["long"] in self.positional_set:
                continue

            option_str = f"  -{opt['short']}, --{opt['long']}" if opt['short'] else f"  --{opt['long']}"
            if opt["has_arg"]:
                arg_help = opt["arg_help"] if opt["arg_help"] else "arg"
                if opt["has_implicit"]:
                    option_str += f" [={arg_help}(={opt['implicit_value']})]"
                else:
                    option_str += f" {arg_help}"

            desc = opt["desc"]
            if opt["has_default"]:
                desc += f" (default: {opt['default_value']})"

            desc_lines = self._wrap_text(desc, longest + 2, 76 - longest - 2)
            result += option_str.ljust(longest + 2) + desc_lines[0] + "\n"
            for line in desc_lines[1:]:
                result += " " * (longest + 2) + line + "\n"

        return result

    def _wrap_text(self, text: str, indent: int, width: int) -> List[str]:
        words = text.split()
        lines = []
        current_line = ""
        for word in words:
            if len(current_line) + 1 + len(word) > width:
                lines.append(current_line)
                current_line = " " * indent + word
            else:
                if current_line:
                    current_line += " " + word
                else:
                    current_line = word
        if current_line:
            lines.append(current_line)
        return lines

    def help(self, groups: List[str] = [""]) -> str:
        result = self.help_string + "\nUsage:\n  " + self.program + " [OPTION...]"
        if self.positional:
            result += " positional parameters"
        result += "\n\n"

        for group in groups:
            group_help = self.help_one_group(group)
            if group_help:
                result += group_help + "\n"

        return result

class OptionAdder:
    def __init__(self, options: Options, group: str):
        self.options = options
        self.group = group

    def __call__(self, opts: str, desc: str, value: Value = StandardValue(False), arg_help: str = "") -> 'OptionAdder':
        match = re.match(r"(([a-zA-Z0-9]),)?([a-zA-Z0-9][-_a-zA-Z0-9]+)", opts)
        if not match:
            raise InvalidOptionFormatError(opts)

        short = match.group(2)
        long = match.group(3)
        self.options.add_option(self.group, short, long, desc, value, arg_help)
        return self

# Example usage
if __name__ == "__main__":
    options = Options("example_program", "Example program help string")
    options.add_options("Main")("o,output", "Output file", StandardValue(str()), "FILE")
    options.add_options("Main")("v,verbose", "Enable verbose logging", StandardValue(bool()))
    options.parse_positional(["output"])

    try:
        options.parse(sys.argv[1:])
        print(options.help(["Main"]))
    except OptionException as e:
        print(e)