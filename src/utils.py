def read_cmd_args(sys_args, commands_str):
    kwargs = dict()
    commands = commands_str.split(" ")
    for i, arg in enumerate(sys_args):
        for command in commands:
            if arg == '-{}'.format(command):
                kwargs[command] = sys_args[i + 1]
    assert len(kwargs) == len(commands), "Command-line args missing, " \
                                         "present=<{}>".format(list(kwargs.keys()))
    return kwargs