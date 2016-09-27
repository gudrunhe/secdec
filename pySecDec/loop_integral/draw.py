'Routines to draw Feynman diagrams'

import subprocess, os

def plot_diagram(internal_lines, external_lines, filename, powerlist=None, neato='neato', extension='pdf', Gstart=0):
    '''
    Draw a Feynman diagram using Graphviz (neato).

    Thanks to Viktor Papara <papara@mpp.mpg.de> for
    his major contributions to this function.

    :param internal_lines:
        list;
        Adjacency list of internal lines,
        e.g. ``[['m',['a',4]],['m',[4,5]],
        ['m',['a',5]],[0,[1,2]],[0,[4,1]],[0,[2,5]]]``

    :param external_lines:
        list;
        Adjacency list of external lines,
        e.g. [['p1',1],['p2',2],['p3','a']]

    :param filename:
        string;
        The name of the output file.
        The generated file gets this name
        plus the file `extension`.

    :param powerlist:
        list, optional;
        The powers of the propagators defined
        by the `internal_lines`.

    :param neato:
        string, default: "neato";
        The shell command to call "neato".

    :param extension:
        string, default: "pdf";
        The file extension. This also defines
        the output format.

    :param Gstart:
        nonnegative int;
        The is value is passed to "neato" with
        the "-Gstart" option. Try changing this
        value if the visualization looks bad.

    '''
    # pinch internal lines
    _pinch(internal_lines, external_lines, powerlist)


    internal_edges_dot = ""
    external_edges_dot = ""
    nodes_dot = ""

    # Write dot language string for internal lines into variable
    # internal_edges_dot.
    dummy = 0
    extra_edges_dummy = [100]
    for edge in internal_lines:
        if powerlist == None:
            power = 1
        else:
            power = powerlist[dummy]
        internal_edges_dot += _edgepower(edge[1][0], edge[1][1], edge[0], power, extra_edges_dummy)
        dummy += 1
        # internal_edges_dot += ( str(edge[1][0])
        #     + " -- "
        #     + str(edge[1][1])
        #     + " [fontname=Arial, labelfloat=true, fontsize=8, label="
        #     + str(edge[0]) + "];\n")

    # Write dot language string for external lines into variable
    # external_edges_dot.
    dummy = 1
    for edge in external_lines:
        external_edges_dot += ( "E" + str(dummy) + " -- "
            + str(edge[1])
            + " [fontname=Arial, labelfloat=true, fontsize=8, taillabel="
            + str(edge[0])
            + "];\n")
        nodes_dot += " E" + str(dummy)
        dummy += 1


    try:
        # Write temporary dot file
        with open(filename, "w") as dot:
            dot.write( "graph 1 {\n" + "node [shape=point];\n"
            "{node [shape=plaintext label=\"\"] " + nodes_dot + "}\n"
            + internal_edges_dot
            + external_edges_dot
            + "}\n"
            )
        # Run neato on temporary dot file
        subprocess.check_call([neato, "-T"+extension, "-Gstart="+str(Gstart), "-Gepsilon=0.000001", "-O", filename])
    finally:
        # Remove temporary file
        os.remove(filename)
        # print internal_lines

def _edgepower(begin, end, label, power=1, uniq=None):
    if power == 1:
        ret = ( str(begin) + " -- " + str(end)
            + " [fontname=Arial, labelfloat=true, fontsize=8, label="
            + str(label) + "];\n")
        return ret
    elif power == 0:
        return ""
    elif power == -1:
        ret = ( str(begin) + " -- " + str(end)
            + " [fontname=Arial, labelfloat=true, fontsize=8, style=dashed, label="
            + str(label) + "];\n")
        return ret
    else:
        # This is for both positive and negative powers

        # Build nodes
        nodes = [str(begin)]
        for i in range(1, abs(power)):
            nodes += [str(uniq[0])]
            uniq[0] += 1
        nodes += [str(end)]

        length = 1./abs(power)
        edges_dot = "edge[len=" + str(length) + "];\n"
        for i in range(0, abs(power)):
            # Call itself to get string for power=+-1
            sign = 1 if ( power > 0) else -1
            edges_dot += _edgepower(nodes[i], nodes[i+1], label, sign)
        edges_dot += "edge[len=1];\n"
        return edges_dot


def _pinch(internal_lines, external_lines, powerlist=None):
    if powerlist != None:
        for idx, power in enumerate(powerlist):
            if power == 0:
                new = internal_lines[idx][1][0]
                old = internal_lines[idx][1][1]
                for line in internal_lines:
                    for i in range(0, len(line[1]) ):
                        if line[1][i] == old:
                            line[1][i] = new
                for line in external_lines:
                    if line[1] == old:
                        line[1] = new
