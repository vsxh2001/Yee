# Print an optspec for argparse to handle cmd's options that are independent of any subcommand.
function __fish_yee_global_optspecs
	string join \n h/help V/version
end

function __fish_yee_needs_command
	# Figure out if the current invocation already has a command.
	set -l cmd (commandline -opc)
	set -e cmd[1]
	argparse -s (__fish_yee_global_optspecs) -- $cmd 2>/dev/null
	or return
	if set -q argv[1]
		# Also print the command, so this can be used to figure out what it is.
		echo $argv[1]
		return 1
	end
	return 0
end

function __fish_yee_using_subcommand
	set -l cmd (__fish_yee_needs_command)
	test -z "$cmd"
	and return 1
	contains -- $cmd[1] $argv
end

complete -c yee -n "__fish_yee_needs_command" -s h -l help -d 'Print help'
complete -c yee -n "__fish_yee_needs_command" -s V -l version -d 'Print version'
complete -c yee -n "__fish_yee_needs_command" -f -a "validate" -d 'Run the validation suite for a given solver (Phase 0: prints planned cases)'
complete -c yee -n "__fish_yee_needs_command" -f -a "mesh" -d 'Mesh a geometry file via Gmsh'
complete -c yee -n "__fish_yee_needs_command" -f -a "run" -d 'Run a simulation defined by a project file (Phase 0 stub)'
complete -c yee -n "__fish_yee_needs_command" -f -a "export" -d 'Export results to Touchstone or HDF5'
complete -c yee -n "__fish_yee_needs_command" -f -a "completions" -d 'Generate a shell completion script on stdout'
complete -c yee -n "__fish_yee_needs_command" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c yee -n "__fish_yee_using_subcommand validate" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c yee -n "__fish_yee_using_subcommand mesh" -s h -l help -d 'Print help'
complete -c yee -n "__fish_yee_using_subcommand run" -s h -l help -d 'Print help'
complete -c yee -n "__fish_yee_using_subcommand export" -l format -d 'Output format' -r -f -a "touchstone\t'Touchstone v1.1 (.s1p/.s2p/.s3p/.s4p)'
hdf5\t'HDF5 (not yet enabled)'"
complete -c yee -n "__fish_yee_using_subcommand export" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c yee -n "__fish_yee_using_subcommand completions" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c yee -n "__fish_yee_using_subcommand help; and not __fish_seen_subcommand_from validate mesh run export completions help" -f -a "validate" -d 'Run the validation suite for a given solver (Phase 0: prints planned cases)'
complete -c yee -n "__fish_yee_using_subcommand help; and not __fish_seen_subcommand_from validate mesh run export completions help" -f -a "mesh" -d 'Mesh a geometry file via Gmsh'
complete -c yee -n "__fish_yee_using_subcommand help; and not __fish_seen_subcommand_from validate mesh run export completions help" -f -a "run" -d 'Run a simulation defined by a project file (Phase 0 stub)'
complete -c yee -n "__fish_yee_using_subcommand help; and not __fish_seen_subcommand_from validate mesh run export completions help" -f -a "export" -d 'Export results to Touchstone or HDF5'
complete -c yee -n "__fish_yee_using_subcommand help; and not __fish_seen_subcommand_from validate mesh run export completions help" -f -a "completions" -d 'Generate a shell completion script on stdout'
complete -c yee -n "__fish_yee_using_subcommand help; and not __fish_seen_subcommand_from validate mesh run export completions help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
