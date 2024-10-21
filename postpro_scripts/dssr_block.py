'''
Run this script with PyMOL, it adds the "dssr_block" command.

(c) Thomas Holder, Schrodinger LLC

License: BSD-2
'''

from pymol import cmd, CmdException

@cmd.extend
def dssr_block(selection='all', state=-1,
        block_file='face',
        block_depth=0.5,
        name='',
        exe='x3dna-dssr',
        quiet=1):
    '''
DESCRIPTION

    Create a nucleid acid cartoon with DSSR

USAGE

    dssr_block [selection [, state [, block_file [, block_depth
        [, name [, exe]]]]]]

ARGUMENTS

    selection = str: atom selection {default: all}

    state = int: object state (0 for all states) {default: -1, current state}

    block_file = face|edge|wc|equal|minor|gray {default: face}

    block_depth = float: thickness of rectangular blocks {default: 0.5}

    name = str: name of new CGO object {default: dssr_block##}

    exe = str: path to "x3dna-dssr" executable {default: x3dna-dssr}

EXAMPLE

    fetch 1ehz, async=0
    as cartoon
    dssr_block
    set cartoon_ladder_radius, 0.1
    set cartoon_ladder_color, gray
    set cartoon_nucleic_acid_mode, 1

    # multi-state
    fetch 2n2d, async=0
    dssr_block 2n2d, 0
    set all_states
    '''
    import subprocess
    import tempfile, os

    state, quiet = int(state), int(quiet)

    tmpfilepdb = tempfile.mktemp('.pdb')
    tmpfiler3d = tempfile.mktemp('.r3d')

    args = [exe,
        '--block-file=' + block_file,
        '--block-depth=' + str(block_depth),
        '-i=' + tmpfilepdb,
        '-o=' + tmpfiler3d,
    ]

    if not name:
        name = cmd.get_unused_name('dssr_block')

    states = [state] if state != 0 else \
            range(1, cmd.count_states(selection) + 1)

    try:
        for state in states:
            cmd.save(tmpfilepdb, selection, state)
            process = subprocess.check_call(args)
            cmd.load(tmpfiler3d, name, max(1, state), zoom=0)
    except subprocess.CalledProcessError:
        raise CmdException('Error: "' + exe + '" failed')
    except OSError:
        raise CmdException('Error: Cannot execute exe="' + exe + '"')
    finally:
        try:
            os.remove(tmpfilepdb)
            os.remove(tmpfiler3d)
        except OSError:
            pass

# tab-completion of arguments
cmd.auto_arg[0].update({
    'dssr_block'    : cmd.auto_arg[0]['zoom'],
})

# vi: expandtab:smarttab
