from yaml import safe_load
import re
import os

import pySecDec as psd

# Helper Function: safely load yaml files
def load_yaml(family_file):
    with open(family_file, 'r') as infile:
        inyaml = safe_load(infile)
    return inyaml

def parse_config(config_dir):

    families_yaml = load_yaml(os.path.join(config_dir,'integralfamilies.yaml'))
    kinematics_yaml = load_yaml(os.path.join(config_dir,'kinematics.yaml'))

    # Extract loop_momenta
    loop_momenta = families_yaml['integralfamilies'][0]['loop_momenta']
    for family in families_yaml['integralfamilies']:
        assert family['loop_momenta'] == loop_momenta, 'Unexpected loop momenta found in integral families' + str(family['loop_momenta'])  + " expected " + str(loop_momenta)

    # Extract independent external momenta
    external_momenta = kinematics_yaml['kinematics']['incoming_momenta'] + kinematics_yaml['kinematics']['outgoing_momenta']
    external_momenta.remove(kinematics_yaml['kinematics']['momentum_conservation'][0]) # remove momenta constrained from momentum conservation

    # Extract replacement_rules
    replacement_rules = kinematics_yaml['kinematics']['scalarproduct_rules']
    replacement_rules = [ ('*'.join(rule[0]), rule[1]) for rule in replacement_rules]

    # Convert propagators to secdec convention
    secdec_integral_families = []
    for integral_family in families_yaml['integralfamilies']:
        secdec_integral_family = {}
        secdec_integral_family['name'] = integral_family['name']
        secdec_propagators = []
        for propagator in integral_family['propagators']:
            # convert ['p1+k1', 'm2'] -> (p1+k1)**2-m2
            if str(propagator[1]) != '0':
                secdec_propagator = '('+str(propagator[0])+')**2-'+str(propagator[1])
            else:
                secdec_propagator = '('+str(propagator[0])+')**2'
            secdec_propagators.append(secdec_propagator)
        secdec_integral_family['propagators'] = secdec_propagators
        secdec_integral_families.append(secdec_integral_family)

    config_dict = {}
    config_dict['loop_momenta'] = loop_momenta
    config_dict['regulators'] = ['eps']
    config_dict['external_momenta'] = external_momenta
    config_dict['replacement_rules'] = replacement_rules
    config_dict['integralfamilies'] = secdec_integral_families
    config_dict['real_parameters'] = [ x for (x,y) in kinematics_yaml['kinematics']['kinematic_invariants']]

    return config_dict

def parse_integral(line, config_dict):
    m = re.match(r' *(\w+)+(( +-?\d+)+)(#.*)?', line)
    fam_name = m.group(1)
    index_list = m.group(2).split(' ')
    index_list = [x for x in index_list if x != '']
    index_list = index_list[4:] # drop tidrs
    # TODO: Compute dimension of integral
    # TODO: Handle crossings
    for family in config_dict['integralfamilies']:
        if family['name'] == fam_name:
            propagators = family['propagators']
    li = psd.loop_integral.LoopIntegralFromPropagators(
            loop_momenta=config_dict['loop_momenta'],
            external_momenta=config_dict['external_momenta'],
            propagators=propagators,
            replacement_rules=config_dict['replacement_rules'],
            powerlist=index_list,
            #dimensionality=str(line_dim) + "-2*eps"
         )
    li_package = psd.LoopPackage(name = fam_name + '_' + '_'.join(index_list).replace('-','m'), loop_integral = li, real_parameters = config_dict['real_parameters'], decomposition_method='geometric')
    #print(line)
    #print(fam_name)
    #print(propagators)
    #print(config_dict['loop_momenta'])
    #print(config_dict['external_momenta'])
    #print(config_dict['replacement_rules'])
    #print(index_list)
    return li_package

def parse_amplitude(amplitude_file, config_dict):

    families = tuple([x['name'] for x in config_dict['integralfamilies']])

    integrals = []
    coefficients = []
    with open(amplitude_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('name') or line == '':
                # Start of amplitude / blank line
                continue
            elif line.startswith(families):
                # Integral
                print('parsing integral:', line)
                li = parse_integral(line, config_dict)
                integrals.append(li)
            elif line.startswith(';'):
                # End of amplitude
                assert len(integrals) == len(coefficients)
                return (integrals, [coefficients])
            else:
                # Coefficient
                print('parsing coefficient')
                coefficients.append(line) 

    return (integrals, [coefficients])

