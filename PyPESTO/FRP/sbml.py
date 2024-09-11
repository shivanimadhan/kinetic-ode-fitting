from libsbml import Model, SBMLDocument, Compartment, Species, Parameter, Reaction, KineticLaw, SpeciesReference, AssignmentRule, ASTNode, parseL3Formula, writeSBMLToString, UNIT_KIND_SECOND, OperationReturnValue_toString, LIBSBML_OPERATION_SUCCESS
import SBMLDiagrams
from typing import Dict, Tuple
import amici
import os

SBML_LEVEL = 3
SBML_VERSION = 2


def check(value, message):
   """If 'value' is None, prints an error message constructed using
   'message' and then exits with status code 1.  If 'value' is an integer,
   it assumes it is a libSBML return status code.  If the code value is
   LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
   prints an error message constructed using 'message' along with text from
   libSBML explaining the meaning of the code, and exits with status code 1.
   """
   if value == None:
     raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
   elif type(value) is int:
     if value == LIBSBML_OPERATION_SUCCESS:
       return
     else:
       err_msg = 'Error encountered trying to ' + message + '.' \
                 + 'LibSBML returned error code ' + str(value) + ': "' \
                 + OperationReturnValue_toString(value).strip() + '"'
       raise SystemExit(err_msg)
   else:
     return
 
def outputSBML(document: SBMLDocument, model_filename: str = 'sbml_model.xml'):
    """Outputs the given model string to the given filename."""
    model_xml_string = writeSBMLToString(document)
    
    with open(model_filename, 'w') as f:
        f.write(model_xml_string)
 
def load_amici_model(model_filepath: str, model_output_dir: str) -> Model:
  
  # Load the model module
  model_name = os.path.basename(model_filepath).split('.')[0]
  model_module = amici.import_model_module(model_name, model_output_dir)

  # Instantiate model
  model = model_module.getModel()

  # Instantiate solver
  solver = model.getSolver()
        
  return model

def load_amici_from_sbml(model_name, load_sbml_model_func):
    
    sbml_model_filepath = f'/SBML/PyPESTO/FRP/{model_name}/sbml_model.xml'
    model_dir = os.path.dirname(sbml_model_filepath)
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    
    sbml_model = load_sbml_model_func(sbml_model_filepath)
    
    print('Importing AMICI model from SBML')
    sbml_importer = amici.SbmlImporter(sbml_model_filepath)
    
    model_output_dir = os.path.join('tmp/', model_name)
    constant_parameters = ["f"]
    observables = {
        "obs_a": {'formula': 'A'},
        "obs_b": {'formula': 'B'},}
    sbml_importer.sbml2amici(
        model_name, 
        model_output_dir, 
        verbose=False,
        observables=observables,
        constant_parameters=constant_parameters
    )
    
    model_module = amici.import_model_module(model_name, model_output_dir)

    # Instantiate model
    model = model_module.getModel()
    
    return model, sbml_model_filepath


def visualize_model(model_filename: str, output_filename: str = 'sbml_model.jpg') -> None:
  
  df = SBMLDiagrams.load(model_filename)

  #df.autolayout()
  df.draw(output_fileName=output_filename)
 
def create_model() -> Tuple[SBMLDocument, Model]:
  
  try:
    document = SBMLDocument(SBML_LEVEL, SBML_VERSION)
  except ValueError:
      raise SystemExit('Could not create SBMLDocumention object')

  model = document.createModel()
  check(model,                              'create model')
  check(model.setTimeUnits("second"),       'set model-wide time units')
  check(model.setExtentUnits("mole"),       'set model units of extent')
  check(model.setSubstanceUnits('mole'),    'set model substance units')
  
  per_second = model.createUnitDefinition()
  check(per_second,                         'create unit definition')
  check(per_second.setId('per_second'),     'set unit definition id')

  unit = per_second.createUnit()
  check(unit,                               'create unit on per_second')
  check(unit.setId('per_second'),           'set unit id')
  check(unit.setKind(UNIT_KIND_SECOND),     'set unit kind')
  check(unit.setExponent(-1),               'set unit exponent')
  check(unit.setScale(0),                   'set unit scale')
  check(unit.setMultiplier(1),              'set unit multiplier')
  
  return document, model

def create_compartment(model: Model, id: str, size: float = 1.0, spatialDimensions: int = 3, units: str = 'litre', isConstant=True) -> Compartment:
  c = model.createCompartment()
  check(c,                                         'create compartment')
  check(c.setId(id),                               'set compartment id')
  check(c.setConstant(isConstant),                 'set compartment "constant"')
  check(c.setSize(size),                           'set compartment "size"')
  check(c.setSpatialDimensions(spatialDimensions), 'set compartment dimensions')
  check(c.setUnits(units),                         'set compartment size units')
  return c
 
def create_species(model: Model, id:str, initialAmount: float = 0.0) -> Species:
  
  s1: Species    = model.createSpecies()
  c: Compartment = model.getCompartment(0)
  
  check(s1,                                 'create species s1')
  check(s1.setId(id),                       'set species s1 id')
  check(s1.setCompartment(c.getId()),            'set species s1 compartment')
  check(s1.setConstant(False),              'set "constant" attribute on s1')
  check(s1.setInitialAmount(initialAmount), 'set initial amount for s1')
  check(s1.setSubstanceUnits('mole'),       'set substance units for s1')
  check(s1.setBoundaryCondition(False),     'set "boundaryCondition" on s1')
  check(s1.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits" on s1')
  return s1
 
def create_parameter(model: Model, id: str, value:float = 0.0, constant: bool = True, units: str = '') -> Parameter:
  k: Parameter = model.createParameter()
  
  # assert UnitKind_isValidUnitKindString(units, level=3, version=2), f'Invalid unit: {units}'
  check(k,                       'create parameter k')
  check(k.setId(id),             'set parameter k id')
  check(k.setConstant(constant), 'set parameter k "constant"')
  check(k.setValue(value),       'set parameter k value')
  check(k.setUnits(units),       'set parameter k units')
  # assert(k.getId() == id)

  return k

def create_reaction(model: Model, id: str, reactantsDict: Dict, productsDict: Dict, kineticLaw: str) -> Reaction:

  r: Reaction = model.createReaction()
  check(r,                                 'create reaction')
  check(r.setId(id),                       'set reaction id')
  check(r.setReversible(False),            'set reaction reversibility flag')
  
  for reactant, stoich in reactantsDict.items():
      species_ref1: SpeciesReference = r.createReactant()
      check(species_ref1,                       'create reactant')
      check(species_ref1.setSpecies(reactant),  'assign reactant species')
      check(species_ref1.setConstant(False),    'set "constant" on species ref 1')
      check(species_ref1.setStoichiometry(stoich), 'set stoichiometry on species ref 1')

  for product, stoich in productsDict.items():
      species_ref2: SpeciesReference = r.createProduct()
      check(species_ref2,                       'create product')
      check(species_ref2.setSpecies(product),   'assign product species')
      check(species_ref2.setConstant(False),     'set "constant" on species ref 2')
      check(species_ref2.setStoichiometry(stoich), 'set stoichiometry on species ref 2')

  # Create kinetic law
  c1: Compartment = model.getCompartment(0)
  
  kineticLaw = f'{kineticLaw} * {c1.getId()}'
  
  math_ast: ASTNode = parseL3Formula(kineticLaw)
  check(math_ast, 'create AST for rate expression')

  kin_law: KineticLaw = r.createKineticLaw()
  check(kin_law,                   'create kinetic law')
  check(kin_law.setMath(math_ast), 'set math on kinetic law')

  return r

def create_rule(model: Model, var, formula: str = '') -> AssignmentRule:
  
  rule: AssignmentRule = model.createAssignmentRule()
  check(rule.setVariable(var.getId()), 'set variable')
  
  # math_ast: ASTNode = parseL3Formula(f'{species1.getId()} + {species2.getId()} + 1e-10')
  math_ast: ASTNode = parseL3Formula(formula)
  check(math_ast,               'create AST for rate expression')
  check(rule.setMath(math_ast), 'set math on kinetic law')
  return rule

def create_SBML_FRP1_v1(model_filepath: str, with_rules: bool = False) -> Model:
  
  print('Creating SBML model (FRP2 v3)')
  
  document, model = create_model()
  c = create_compartment(model, 'c')
  
  print('Creating species.')
  I   = create_species(model, 'I', initialAmount=0.005)
  R   = create_species(model, 'R')

  M   = create_species(model, 'M', initialAmount=1.5)
  P   = create_species(model, 'P')
  PD  = create_species(model, 'PD')

  print('Generating parameters.')
  
  kp = create_parameter(model, 'kp', value=1, constant=False)
  kp_kt_ratio = create_parameter(model, 'kp_kt_ratio', value=1, constant=False)
  kd_kt = create_parameter(model, 'kd_kt', value=1, constant=False)
  Kr = create_parameter(model, 'Kr', value=1, constant=False)
  # kp_kt_ratio_true = kp_true / kt_true
# kd_kt_true = kd_true * kt_true
  
  kd = create_parameter(model, 'kd', value=3e-6, constant=False)
  f  = create_parameter(model, 'f', constant=True, value=0.5)
  # kp = create_parameter(model, 'kp', value=1, constant=False)
  kt = create_parameter(model, 'kt', value=0, constant=False)
  kdp = create_parameter(model, 'kdp', value=1, constant=False)
  
  create_rule(model, kt, formula=f'{kp.getId()} / {kp_kt_ratio.getId()}')
  create_rule(model, kd, formula=f'{kd_kt.getId()} / {kt.getId()}')
  create_rule(model, kdp, formula=f'{Kr.getId()} * {kp.getId()}')

  
  # Defining reactions
  # Syntax: (reaction_id, {reactants: stoich}, {products: stoich}, kinetic_law)
  reactions = [
      ('initiation', {I.getId(): 1}, {R.getId(): 2*f.getValue()}, f'{kd.getId()} * {I.getId()}'),
      ('init_prop', {R.getId(): 1, M.getId(): 1}, {P.getId(): 1}, f'{kp.getId()} * {R.getId()} * {M.getId()}'),
      ('prop', {P.getId(): 1, M.getId(): 1}, {P.getId(): 1}, f'{kp.getId()} * {P.getId()} * {M.getId()}'),
      ('deprop', {P.getId(): 1}, {P.getId(): 1, M.getId(): 1}, f'{kdp.getId()} * {P.getId()}'),
      ('term', {P.getId(): 2}, {PD.getId(): 1}, f'{kt.getId()} * {P.getId()} * {P.getId()}'),
  ]
  # create_rule(model, ktAB, formula=f'({ktAA.getId()} * {ktBB.getId()})^0.5')
  
  # (reaction_id, reactants_dict, products_dict, kinetic_law)
  print('Creating reactions')
  generated_reactions = []
  for r in reactions:
      reaction = create_reaction(model, r[0], r[1], r[2], r[3])
      generated_reactions.append(reaction)
      
  outputSBML(document, model_filepath)
  
  return model


def create_SBML_FRP2_v1(model_filepath: str, with_rules: bool = False) -> Model:
  
  print('Creating SBML model (FRP2 v1)')
  
  document, model = create_model()
  c = create_compartment(model, 'c')
  
  print('Creating species.')
  I   = create_species(model, 'I', initialAmount=0.000)
  R   = create_species(model, 'R', initialAmount=0.005)

  A   = create_species(model, 'A', initialAmount=1.5)
  B   = create_species(model, 'B', initialAmount=1.5)

  RA  = create_species(model, 'RA')
  RB  = create_species(model, 'RB')

  PAA = create_species(model, 'PAA')
  PAB = create_species(model, 'PAB')
  PBA = create_species(model, 'PBA')
  PBB = create_species(model, 'PBB')

  print('Generating parameters.')
  kd = create_parameter(model, 'kd', value=1e10)
  f  = create_parameter(model, 'f', constant=True, value=0.5)
  kpAA = create_parameter(model, 'kpAA', value=1, constant=False)
  kpAB = create_parameter(model, 'kpAB', value=1, constant=False)
  kpBA = create_parameter(model, 'kpBA', value=1, constant=False)
  kpBB = create_parameter(model, 'kpBB', value=1, constant=False)
  
  rA = create_parameter(model, 'rA', value=1, constant=False) # rA = kpAA / kpAB
  rB = create_parameter(model, 'rB', value=1, constant=False) # rB = kpBB / kpBA
  rX = create_parameter(model, 'rX', value=1, constant=False) # rX = kpAA / kpBB
  
  # Defining reactions
  # Syntax: (reaction_id, {reactants: stoich}, {products: stoich}, kinetic_law)
  reactions = [
      ('initiation', {I.getId(): 1}, {R.getId(): 2*f.getValue()}, f'{kd.getId()} * {I.getId()}'),
      ('prop_A', {R.getId(): 1, A.getId(): 1}, {RA.getId(): 1}, f'{kpAA.getId()} * {R.getId()} * {A.getId()}'),
      ('prop_B', {R.getId(): 1, B.getId(): 1}, {RB.getId(): 1}, f'{kpBB.getId()} * {R.getId()} * {B.getId()}'),
      ('prop_RA_A', {RA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {RA.getId()} * {A.getId()}'),
      ('prop_RA_B', {RA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {RA.getId()} * {B.getId()}'),
      ('prop_RB_A', {RB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {RB.getId()} * {A.getId()}'),
      ('prop_RB_B', {RB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {RB.getId()} * {B.getId()}'),
      ('prop_PAA_A', {PAA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {PAA.getId()} * {A.getId()}'),
      ('prop_PAA_B', {PAA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {PAA.getId()} * {B.getId()}'),
      ('prop_PAB_A', {PAB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {PAB.getId()} * {A.getId()}'),
      ('prop_PAB_B', {PAB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {PAB.getId()} * {B.getId()}'),
      ('prop_PBA_A', {PBA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {PBA.getId()} * {A.getId()}'),
      ('prop_PBA_B', {PBA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {PBA.getId()} * {B.getId()}'),
      ('prop_PBB_A', {PBB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {PBB.getId()} * {A.getId()}'),
      ('prop_PBB_B', {PBB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {PBB.getId()} * {B.getId()}'),
  ]

  # Define assignment rules
  create_rule(model, kpAB, formula=f'{kpAA.getId()} / {rA.getId()}') # rA = kpAA / kpAB -> kpAB = kpAA / rA
  create_rule(model, kpBB, formula=f'{kpAA.getId()} / {rX.getId()}') # rX = kpAA / kpBB -> kpBB = kpAA / rX
  create_rule(model, kpBA, formula=f'{kpBB.getId()} / {rB.getId()}') # rB = kpBB / kpBA -> kpBA = kpBB / rB

  # (reaction_id, reactants_dict, products_dict, kinetic_law)
  print('Creating reactions')
  generated_reactions = []
  for r in reactions:
      reaction = create_reaction(model, r[0], r[1], r[2], r[3])
      generated_reactions.append(reaction)
      
  outputSBML(document, model_filepath)
  
  return model

def create_SBML_FRP2_v2(model_filepath: str, with_rules: bool = False) -> Model:
  
  print('Creating SBML model (FRP2 v2)')
  
  document, model = create_model()
  c = create_compartment(model, 'c')
  
  print('Creating species.')
  I   = create_species(model, 'I', initialAmount=0.000)
  R   = create_species(model, 'R', initialAmount=0.005)

  A   = create_species(model, 'A', initialAmount=1.5)
  B   = create_species(model, 'B', initialAmount=1.5)

  RA  = create_species(model, 'RA')
  RB  = create_species(model, 'RB')

  PAA = create_species(model, 'PAA')
  PAB = create_species(model, 'PAB')
  PBA = create_species(model, 'PBA')
  PBB = create_species(model, 'PBB')
  
  PA = create_species(model, 'PA')
  PB = create_species(model, 'PB')

  print('Generating parameters.')
  kd = create_parameter(model, 'kd', value=1e10)
  f  = create_parameter(model, 'f', constant=True, value=0.5)
  kpAA = create_parameter(model, 'kpAA', value=1, constant=False)
  kpAB = create_parameter(model, 'kpAB', value=1, constant=False)
  kpBA = create_parameter(model, 'kpBA', value=1, constant=False)
  kpBB = create_parameter(model, 'kpBB', value=1, constant=False)
  kdAA = create_parameter(model, 'kdAA', value=1, constant=False)
  
  rA = create_parameter(model, 'rA', value=1, constant=False) # rA = kpAA / kpAB
  rB = create_parameter(model, 'rB', value=1, constant=False) # rB = kpBB / kpBA
  rX = create_parameter(model, 'rX', value=1, constant=False) # rX = kpAA / kpBB
  
  KAA = create_parameter(model, 'KAA', value=0, constant=False) # KAA = kpAA / kdAA
  # KAB = create_parameter(model, 'rAB', value=0, constant=False) # KAB = kpAB / kdAB
  # KBA = create_parameter(model, 'rBA', value=0, constant=False) # KBA = kpBA / kdBA
  # KBB = create_parameter(model, 'rBB', value=0, constant=False) # KBB = kpBB / kdBB
  
  fPAA = create_parameter(model, 'fPAA', value=1, constant=False) # fPAA = PAA / PA
  fPAB = create_parameter(model, 'fPAB', value=1, constant=False) # fPAB = PAB / PA
  fPBA = create_parameter(model, 'fPBA', value=1, constant=False) # fPBA = PBA / PB
  fPBB = create_parameter(model, 'fPBB', value=1, constant=False) # fPBB = PBB / PB
  
  create_rule(model, PA, formula=f'{PAA.getId()} + {PBA.getId()}')
  create_rule(model, PB, formula=f'{PAB.getId()} + {PBB.getId()}')
  
  eps = 1e-10
  create_rule(model, fPAA, formula=f'({PAA.getId()} + {eps}) / ({PA.getId()} + {eps})')
  create_rule(model, fPAB, formula=f'({PAB.getId()} + {eps}) / ({PB.getId()} + {eps})')
  create_rule(model, fPBA, formula=f'({PBA.getId()} + {eps}) / ({PA.getId()} + {eps})')
  create_rule(model, fPBB, formula=f'({PBB.getId()} + {eps}) / ({PB.getId()} + {eps})')
  
  # Defining reactions
  # Syntax: (reaction_id, {reactants: stoich}, {products: stoich}, kinetic_law)
  reactions = [
      ('initiation', {I.getId(): 1}, {R.getId(): 2*f.getValue()}, f'{kd.getId()} * {I.getId()}'),
      ('prop_A', {R.getId(): 1, A.getId(): 1}, {RA.getId(): 1}, f'{kpAA.getId()} * {R.getId()} * {A.getId()}'),
      ('prop_B', {R.getId(): 1, B.getId(): 1}, {RB.getId(): 1}, f'{kpBB.getId()} * {R.getId()} * {B.getId()}'),
      ('prop_RA_A', {RA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {RA.getId()} * {A.getId()}'),
      ('prop_RA_B', {RA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {RA.getId()} * {B.getId()}'),
      ('prop_RB_A', {RB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {RB.getId()} * {A.getId()}'),
      ('prop_RB_B', {RB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {RB.getId()} * {B.getId()}'),
      ('prop_PAA_A', {PAA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {PAA.getId()} * {A.getId()}'),
      ('prop_PAA_B', {PAA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {PAA.getId()} * {B.getId()}'),
      ('prop_PAB_A', {PAB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {PAB.getId()} * {A.getId()}'),
      ('prop_PAB_B', {PAB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {PAB.getId()} * {B.getId()}'),
      ('prop_PBA_A', {PBA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {PBA.getId()} * {A.getId()}'),
      ('prop_PBA_B', {PBA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {PBA.getId()} * {B.getId()}'),
      ('prop_PBB_A', {PBB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {PBB.getId()} * {A.getId()}'),
      ('prop_PBB_B', {PBB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {PBB.getId()} * {B.getId()}'),
      ('deprop_PAAA', {PAA.getId(): 1}, {PAA.getId(): 1, A.getId(): 1}, f'{kdAA.getId()} * {fPAA.getId()} * {PAA.getId()}'),
      ('deprop_PBAA', {PAA.getId(): 1}, {PBA.getId(): 1, A.getId(): 1}, f'{kdAA.getId()} * {fPBA.getId()} * {PAA.getId()}'),
      # ('deprop_PABA', {PBA.getId(): 1}, {PAB.getId(): 1, A.getId(): 1}, f'{kdBA.getId()} * {fPAB.getId()} * {PBA.getId()}'),
      # ('deprop_PBBA', {PBA.getId(): 1}, {PBB.getId(): 1, A.getId(): 1}, f'{kdBA.getId()} * {fPBB.getId()} * {PBA.getId()}'),
      # ('deprop_PAAB', {PAB.getId(): 1}, {PAA.getId(): 1, B.getId(): 1}, f'{kdAB.getId()} * {fPAA.getId()} * {PAB.getId()}'),
      # ('deprop_PBAB', {PAB.getId(): 1}, {PBA.getId(): 1, B.getId(): 1}, f'{kdAB.getId()} * {fPBA.getId()} * {PAB.getId()}'),
      # ('deprop_PABB', {PBB.getId(): 1}, {PAB.getId(): 1, B.getId(): 1}, f'{kdBB.getId()} * {fPAB.getId()} * {PBB.getId()}'),
      # ('deprop_PBBB', {PBB.getId(): 1}, {PBB.getId(): 1, B.getId(): 1}, f'{kdBB.getId()} * {fPBB.getId()} * {PBB.getId()}'),
  ]
  # assert()
  
  
  # Define assignment rules
  create_rule(model, kpAB, formula=f'{kpAA.getId()} / {rA.getId()}') # rA = kpAA / kpAB -> kpAB = kpAA / rA
  create_rule(model, kpBB, formula=f'{kpAA.getId()} / {rX.getId()}') # rX = kpAA / kpBB -> kpBB = kpAA / rX
  create_rule(model, kpBA, formula=f'{kpBB.getId()} / {rB.getId()}') # rB = kpBB / kpBA -> kpBA = kpBB / rB
  create_rule(model, kdAA, formula=f'{kpAA.getId()} / {KAA.getId()}') # KAA = kpAA / kdAA -> kdAA = kpAA / KAA
  # create_rule(model, kdAB, formula=f'{kpAB.getId()} / {KAB.getId()}') # KAB = kpAB / kdAB -> kdAB = kpAB / KAB
  # create_rule(model, kdBA, formula=f'{kpBA.getId()} / {KBA.getId()}') # KBA = kpBA / kdBA -> kdBA = kpBA / KBA
  # create_rule(model, kdBB, formula=f'{kpBB.getId()} / {KBB.getId()}') # KBB = kpBB / kdBB -> kdBB = kpBB / KBB
  
  # (reaction_id, reactants_dict, products_dict, kinetic_law)
  print('Creating reactions')
  generated_reactions = []
  for r in reactions:
      reaction = create_reaction(model, r[0], r[1], r[2], r[3])
      generated_reactions.append(reaction)
      
  outputSBML(document, model_filepath)
  
  return model

def create_SBML_FRP2_v3(model_filepath: str, with_rules: bool = False) -> Model:
  
  print('Creating SBML model (FRP2 v3)')
  
  document, model = create_model()
  c = create_compartment(model, 'c')
  
  print('Creating species.')
  I   = create_species(model, 'I', initialAmount=0.005)
  R   = create_species(model, 'R', initialAmount=0.000)

  A   = create_species(model, 'A', initialAmount=1.5)
  B   = create_species(model, 'B', initialAmount=1.5)

  RA  = create_species(model, 'RA')
  RB  = create_species(model, 'RB')

  PAA = create_species(model, 'PAA')
  PAB = create_species(model, 'PAB')
  PBA = create_species(model, 'PBA')
  PBB = create_species(model, 'PBB')
  
  PA = create_species(model, 'PA')
  PB = create_species(model, 'PB')
  
  PD = create_species(model, 'PD')

  print('Generating parameters.')
  kpAA = create_parameter(model, 'kpAA', value=1, constant=False)
  kpAB = create_parameter(model, 'kpAB', value=1, constant=False)
  kpBA = create_parameter(model, 'kpBA', value=1, constant=False)
  kpBB = create_parameter(model, 'kpBB', value=1, constant=False)
  kdBB = create_parameter(model, 'kdBB', value=1, constant=False)
  
  rA = create_parameter(model, 'rA', value=1, constant=False) # rA = kpAA / kpAB
  rB = create_parameter(model, 'rB', value=1, constant=False) # rB = kpBB / kpBA
  rX = create_parameter(model, 'rX', value=1, constant=False) # rX = kpAA / kpBB
  KBB = create_parameter(model, 'KBB', value=0, constant=False) # KBB = kpBB / kdBB
  
  kt = create_parameter(model, 'kt', value=0, constant=False)
  kd = create_parameter(model, 'kd', value=3e-6, constant=False)
  f  = create_parameter(model, 'f', constant=True, value=0.5)
  kp_kt_ratio = create_parameter(model, 'kp_kt_ratio', value=1, constant=False)
  kd_kt = create_parameter(model, 'kd_kt', value=1, constant=False)
  
  create_rule(model, kt, formula=f'{kpAA.getId()} / {kp_kt_ratio.getId()}')
  create_rule(model, kd, formula=f'{kd_kt.getId()} / {kt.getId()}')

  # Define assignment rules
  create_rule(model, kpAB, formula=f'{kpAA.getId()} / {rA.getId()}') # rA = kpAA / kpAB -> kpAB = kpAA / rA
  create_rule(model, kpBB, formula=f'{kpAA.getId()} / {rX.getId()}') # rX = kpAA / kpBB -> kpBB = kpAA / rX
  create_rule(model, kpBA, formula=f'{kpBB.getId()} / {rB.getId()}') # rB = kpBB / kpBA -> kpBA = kpBB / rB
  create_rule(model, kdBB, formula=f'{kpBB.getId()} * {KBB.getId()}') # KBB = kpBB / kdBB -> kdBB = kpBB / KBB
  
  
  # kt = create_parameter(model, 'kt', value=0, constant=False)
  # ktAA = create_parameter(model, 'ktAA', value=0, constant=False)
  # ktBB = create_parameter(model, 'ktBB', value=0, constant=False)
  # ktAB = create_parameter(model, 'ktAB', value=0, constant=False)
  
  fPAA = create_parameter(model, 'fPAA', value=1, constant=False) # fPAA = PAA / PA
  fPAB = create_parameter(model, 'fPAB', value=1, constant=False) # fPAB = PAB / PA
  fPBA = create_parameter(model, 'fPBA', value=1, constant=False) # fPBA = PBA / PB
  fPBB = create_parameter(model, 'fPBB', value=1, constant=False) # fPBB = PBB / PB
  
  create_rule(model, PA, formula=f'{PAA.getId()} + {PBA.getId()}')
  create_rule(model, PB, formula=f'{PAB.getId()} + {PBB.getId()}')
  
  eps = 1e-10
  create_rule(model, fPAA, formula=f'({PAA.getId()} + {eps}) / ({PA.getId()} + {eps})')
  create_rule(model, fPAB, formula=f'({PAB.getId()} + {eps}) / ({PB.getId()} + {eps})')
  create_rule(model, fPBA, formula=f'({PBA.getId()} + {eps}) / ({PA.getId()} + {eps})')
  create_rule(model, fPBB, formula=f'({PBB.getId()} + {eps}) / ({PB.getId()} + {eps})')
  
  # Defining reactions
  # Syntax: (reaction_id, {reactants: stoich}, {products: stoich}, kinetic_law)
  reactions = [
      ('initiation', {I.getId(): 1}, {R.getId(): 2*f.getValue()}, f'{kd.getId()} * {f.getId()} * {I.getId()}'),
      ('prop_A', {R.getId(): 1, A.getId(): 1}, {RA.getId(): 1}, f'{kpAA.getId()} * {R.getId()} * {A.getId()}'),
      ('prop_B', {R.getId(): 1, B.getId(): 1}, {RB.getId(): 1}, f'{kpBB.getId()} * {R.getId()} * {B.getId()}'),
      ('prop_RA_A', {RA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {RA.getId()} * {A.getId()}'),
      ('prop_RA_B', {RA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {RA.getId()} * {B.getId()}'),
      ('prop_RB_A', {RB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {RB.getId()} * {A.getId()}'),
      ('prop_RB_B', {RB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {RB.getId()} * {B.getId()}'),
      ('prop_PAA_A', {PAA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {PAA.getId()} * {A.getId()}'),
      ('prop_PAA_B', {PAA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {PAA.getId()} * {B.getId()}'),
      ('prop_PAB_A', {PAB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {PAB.getId()} * {A.getId()}'),
      ('prop_PAB_B', {PAB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {PAB.getId()} * {B.getId()}'),
      ('prop_PBA_A', {PBA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {PBA.getId()} * {A.getId()}'),
      ('prop_PBA_B', {PBA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {PBA.getId()} * {B.getId()}'),
      ('prop_PBB_A', {PBB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {PBB.getId()} * {A.getId()}'),
      ('prop_PBB_B', {PBB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {PBB.getId()} * {B.getId()}'),
      # ('deprop_PAAA', {PAA.getId(): 1}, {PAA.getId(): 1, A.getId(): 1}, f'{kdAA.getId()} * {fPAA.getId()} * {PAA.getId()}'),
      # ('deprop_PBAA', {PAA.getId(): 1}, {PBA.getId(): 1, A.getId(): 1}, f'{kdAA.getId()} * {fPBA.getId()} * {PAA.getId()}'),
      ('deprop_PABB', {PBB.getId(): 1}, {PAB.getId(): 1, B.getId(): 1}, f'{kdBB.getId()} * {fPAB.getId()} * {PBB.getId()}'),
      ('deprop_PBBB', {PBB.getId(): 1}, {PBB.getId(): 1, B.getId(): 1}, f'{kdBB.getId()} * {fPBB.getId()} * {PBB.getId()}'),
      # ('termc_PA_PA', {PA.getId(): 1, PA.getId(): 1}, {PD.getId(): 1}, f'{ktcAA.getId()} * {PA.getId()} * {PA.getId()}'),
      # ('termc_PA_PB', {PA.getId(): 1, PB.getId(): 1}, {PD.getId(): 1}, f'{ktcAB.getId()} * {PA.getId()} * {PB.getId()}'),
      # ('termc_PB_PB', {PB.getId(): 1, PB.getId(): 1}, {PD.getId(): 1}, f'{ktcBB.getId()} * {PB.getId()} * {PB.getId()}'),
      # ('termd_PA_PA', {PA.getId(): 1, PA.getId(): 1}, {PD.getId(): 1}, f'{ktdAA.getId()} * {PA.getId()} * {PA.getId()}'),
      # ('termd_PA_PB', {PA.getId(): 1, PB.getId(): 1}, {PD.getId(): 1}, f'{ktdAB.getId()} * {PA.getId()} * {PB.getId()}'),
      # ('termd_PB_PB', {PB.getId(): 1, PB.getId(): 1}, {PD.getId(): 1}, f'{ktdBB.getId()} * {PB.getId()} * {PB.getId()}'),
      
      # ('term_PAA_PAA', {PAA.getId(): 2}, {PD.getId(): 1}, f'{kt.getId()} * {PAA.getId()} * {PAA.getId()}'),
      # ('term_PAA_PBA', {PAA.getId(): 1, PBA.getId(): 1}, {PD.getId(): 1}, f'{kt.getId()} * {PAA.getId()} * {PBA.getId()}'),
      # ('term_PAA_PAB', {PAA.getId(): 1, PAB.getId(): 1}, {PD.getId(): 1}, f'{kt.getId()} * {PAA.getId()} * {PAB.getId()}'),
      
      
      
      # ('term_PBA_PBA', {PBA.getId(): 2}, {PD.getId(): 1}, f'{kt.getId()} * {PBA.getId()} * {PBA.getId()}'),
      
      # ('term_PAA_PAB', {PAA.getId(): 1, PAB.getId(): 1}, {PD.getId(): 1}, f'{kt.getId()} * {PAA.getId()} * {PAB.getId()}'),
      # ('term_PAA_PBB', {PAA.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{kt.getId()} * {PAA.getId()} * {PBB.getId()}'),
      # ('term_PBA_PAB', {PBA.getId(): 1, PAB.getId(): 1}, {PD.getId(): 1}, f'{kt.getId()} * {PBA.getId()} * {PAB.getId()}'),
      # ('term_PBA_PBB', {PBA.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{kt.getId()} * {PBA.getId()} * {PBB.getId()}'),
      
      # ('term_PAB_PAB', {PAB.getId(): 2}, {PD.getId(): 1}, f'{kt.getId()} * {PAB.getId()} * {PAB.getId()}'),
      # ('term_PAB_PBB', {PAB.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{kt.getId()} * {PAB.getId()} * {PBB.getId()}'),
      # ('term_PBB_PBB', {PBB.getId(): 2}, {PD.getId(): 1}, f'{kt.getId()} * {PBB.getId()} * {PBB.getId()}'),
      
      
      # ('termc_PAA_PAA', {PAA.getId(): 1, PAA.getId(): 1}, {PD.getId(): 1}, f'{ktcAA.getId()} * {PAA.getId()} * {PAA.getId()}'),
      # ('termc_PAA_PBA', {PAA.getId(): 1, PBA.getId(): 1}, {PD.getId(): 1}, f'{ktcAA.getId()} * {PAA.getId()} * {PBA.getId()}'),
      # ('termc_PBA_PBA', {PBA.getId(): 1, PBA.getId(): 1}, {PD.getId(): 1}, f'{ktcAA.getId()} * {PBA.getId()} * {PBA.getId()}'),
      
      # ('termc_PAA_PAB', {PAA.getId(): 1, PAB.getId(): 1}, {PD.getId(): 1}, f'{ktcAB.getId()} * {PAA.getId()} * {PAB.getId()}'),
      # ('termc_PAA_PBB', {PAA.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{ktcAB.getId()} * {PAA.getId()} * {PBB.getId()}'),
      # ('termc_PBA_PAB', {PBA.getId(): 1, PAB.getId(): 1}, {PD.getId(): 1}, f'{ktcAB.getId()} * {PBA.getId()} * {PAB.getId()}'),
      # ('termc_PBA_PBB', {PBA.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{ktcAB.getId()} * {PBA.getId()} * {PBB.getId()}'),
      
      # ('termc_PAB_PAB', {PAB.getId(): 1, PAB.getId(): 1}, {PD.getId(): 1}, f'{ktcBB.getId()} * {PAB.getId()} * {PAB.getId()}'),
      # ('termc_PAB_PBB', {PAB.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{ktcBB.getId()} * {PAB.getId()} * {PBB.getId()}'),
      # ('termc_PBB_PBB', {PBB.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{ktcBB.getId()} * {PBB.getId()} * {PBB.getId()}'),
      
      # ('termd_PAA_PAA', {PAA.getId(): 1, PAA.getId(): 1}, {PD.getId(): 1}, f'{ktdAA.getId()} * {PAA.getId()} * {PAA.getId()}'),
      # ('termd_PAA_PBA', {PAA.getId(): 1, PBA.getId(): 1}, {PD.getId(): 1}, f'{ktdAA.getId()} * {PAA.getId()} * {PBA.getId()}'),
      # ('termd_PBA_PBA', {PBA.getId(): 1, PBA.getId(): 1}, {PD.getId(): 1}, f'{ktdAA.getId()} * {PBA.getId()} * {PBA.getId()}'),
      
      # ('termd_PAA_PAB', {PAA.getId(): 1, PAB.getId(): 1}, {PD.getId(): 1}, f'{ktdAB.getId()} * {PAA.getId()} * {PAB.getId()}'),
      # ('termd_PAA_PBB', {PAA.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{ktdAB.getId()} * {PAA.getId()} * {PBB.getId()}'),
      # ('termd_PBA_PAB', {PBA.getId(): 1, PAB.getId(): 1}, {PD.getId(): 1}, f'{ktdAB.getId()} * {PBA.getId()} * {PAB.getId()}'),
      
      # ('termd_PAB_PAB', {PAB.getId(): 1, PAB.getId(): 1}, {PD.getId(): 1}, f'{ktdBB.getId()} * {PAB.getId()} * {PAB.getId()}'),
      # ('termd_PAB_PBB', {PAB.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{ktdBB.getId()} * {PAB.getId()} * {PBB.getId()}'),
      # ('termd_PBB_PBB', {PBB.getId(): 1, PBB.getId(): 1}, {PD.getId(): 1}, f'{ktdBB.getId()} * {PBB.getId()} * {PBB.getId()}'),
      
      
      
  ]
  from itertools import combinations_with_replacement, permutations
  termination_species = [PAA, PAB, PBA, PBB]
  x = 0
  for s1, s2 in combinations_with_replacement(termination_species, 2):
    reaction_id = f'term_{s1.getId()}_{s2.getId()}'
    reactants = {s1.getId(): 1, s2.getId(): 1}
    products = {PD.getId(): 1}
    
    if (s1.getId() != s2.getId()):
      kinetic_law = f'4 * {kt.getId()} * {s1.getId()} * {s2.getId()}'
    else:
      kinetic_law = f'2 * {kt.getId()} * {s1.getId()} * {s2.getId()}'
      
    reactions.append((reaction_id, reactants, products, kinetic_law))
    x += 1
    
    # if (s1.getId() != s2.getId()):
    #   reaction_id = f'term_{s2.getId()}_{s1.getId()}'
    #   reactants = {s2.getId(): 1, s1.getId(): 1}
    #   products = {PD.getId(): 1}
    #   kinetic_law = f'2 * {kt.getId()} * {s2.getId()} * {s1.getId()}'
    #   reactions.append((reaction_id, reactants, products, kinetic_law))
    #   x += 1
      
  print(f'Num Termination Reactions: {x}')
  
  
  # create_rule(model, ktAB, formula=f'({ktAA.getId()} * {ktBB.getId()})^0.5')
  
  # (reaction_id, reactants_dict, products_dict, kinetic_law)
  print('Creating reactions')
  generated_reactions = []
  for r in reactions:
      reaction = create_reaction(model, r[0], r[1], r[2], r[3])
      generated_reactions.append(reaction)
      
  outputSBML(document, model_filepath)
  
  return model


def create_SBML_FRP3_v1(model_filepath: str, with_rules: bool = False) -> Model:
  
  print('Creating SBML model (FRP3 v1)')
  
  document, model = create_model()
  c = create_compartment(model, 'c')
  
  print('Creating species.')
  I   = create_species(model, 'I', initialAmount=0.000)
  R   = create_species(model, 'R', initialAmount=0.005)

  A   = create_species(model, 'A', initialAmount=1.5)
  B   = create_species(model, 'B', initialAmount=1.5)
  C   = create_species(model, 'C', initialAmount=1.5)

  RA  = create_species(model, 'RA')
  RB  = create_species(model, 'RB')
  RC  = create_species(model, 'RC')

  PAA = create_species(model, 'PAA')
  PAB = create_species(model, 'PAB')
  PAC = create_species(model, 'PAC')
  
  PBA = create_species(model, 'PBA')
  PBB = create_species(model, 'PBB')
  PBC = create_species(model, 'PBC')
  
  PCA = create_species(model, 'PCA')
  PCB = create_species(model, 'PCB')
  PCC = create_species(model, 'PCC')

  print('Generating parameters.')
  kd = create_parameter(model, 'kd', value=1e10)
  f  = create_parameter(model, 'f', constant=True, value=0.5)
  
  kpAA = create_parameter(model, 'kpAA', value=1, constant=False)
  kpAB = create_parameter(model, 'kpAB', value=1, constant=False)
  kpAC = create_parameter(model, 'kpAC', value=1, constant=False)
  
  kpBA = create_parameter(model, 'kpBA', value=1, constant=False)
  kpBB = create_parameter(model, 'kpBB', value=1, constant=False)
  kpBC = create_parameter(model, 'kpBC', value=1, constant=False)
  
  kpCA = create_parameter(model, 'kpCA', value=1, constant=False)
  kpCB = create_parameter(model, 'kpCB', value=1, constant=False)
  kpCC = create_parameter(model, 'kpCC', value=1, constant=False)
  
  rAB = create_parameter(model, 'rAB', value=1, constant=False) # rAB = kpAA / kpAB
  rAC = create_parameter(model, 'rAC', value=1, constant=False) # rAC = kpAA / kpAC
  
  rBA = create_parameter(model, 'rBA', value=1, constant=False) # rBA = kpBB / kpBA
  rBC = create_parameter(model, 'rBC', value=1, constant=False) # rBC = kpBB / kpBC
  
  rCA = create_parameter(model, 'rCA', value=1, constant=False) # rCA = kpCC / kpCA
  rCB = create_parameter(model, 'rCB', value=1, constant=False) # rCB = kpCC / kpCB 
  
  rAxB = create_parameter(model, 'rAxB', value=1, constant=False) # rAxB = kpAA / kpBB
  rAxC = create_parameter(model, 'rAxC', value=1, constant=False) # rAxC = kpAA / kpCC
  
  # Defining reactions
  # Syntax: (reaction_id, {reactants: stoich}, {products: stoich}, kinetic_law)
  reactions = [
      ('initiation', {I.getId(): 1}, {R.getId(): 2}, f'{kd.getId()} * {I.getId()}'),
      ('prop_A', {R.getId(): 1, A.getId(): 1}, {RA.getId(): 1}, f'{kpAA.getId()} * {R.getId()} * {A.getId()}'),
      ('prop_B', {R.getId(): 1, B.getId(): 1}, {RB.getId(): 1}, f'{kpBB.getId()} * {R.getId()} * {B.getId()}'),
      ('prop_C', {R.getId(): 1, C.getId(): 1}, {RC.getId(): 1}, f'{kpCC.getId()} * {R.getId()} * {C.getId()}'),
      ('prop_RA_A', {RA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {RA.getId()} * {A.getId()}'),
      ('prop_RA_B', {RA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {RA.getId()} * {B.getId()}'),
      ('prop_RA_C', {RA.getId(): 1, C.getId(): 1}, {PAC.getId(): 1}, f'{kpAC.getId()} * {RA.getId()} * {C.getId()}'),
      ('prop_RB_A', {RB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {RB.getId()} * {A.getId()}'),
      ('prop_RB_B', {RB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {RB.getId()} * {B.getId()}'),
      ('prop_RB_C', {RB.getId(): 1, C.getId(): 1}, {PBC.getId(): 1}, f'{kpBC.getId()} * {RB.getId()} * {C.getId()}'),
      ('prop_RC_A', {RC.getId(): 1, A.getId(): 1}, {PCA.getId(): 1}, f'{kpCA.getId()} * {RC.getId()} * {A.getId()}'),
      ('prop_RC_B', {RC.getId(): 1, B.getId(): 1}, {PCB.getId(): 1}, f'{kpCB.getId()} * {RC.getId()} * {B.getId()}'),
      ('prop_RC_C', {RC.getId(): 1, C.getId(): 1}, {PCC.getId(): 1}, f'{kpCC.getId()} * {RC.getId()} * {C.getId()}'),
      ('prop_PAA_A', {PAA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {PAA.getId()} * {A.getId()}'),
      ('prop_PAA_B', {PAA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {PAA.getId()} * {B.getId()}'),
      ('prop_PAA_C', {PAA.getId(): 1, C.getId(): 1}, {PAC.getId(): 1}, f'{kpAC.getId()} * {PAA.getId()} * {C.getId()}'),
      ('prop_PAB_A', {PAB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {PAB.getId()} * {A.getId()}'),
      ('prop_PAB_B', {PAB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {PAB.getId()} * {B.getId()}'),
      ('prop_PAB_C', {PAB.getId(): 1, C.getId(): 1}, {PBC.getId(): 1}, f'{kpBC.getId()} * {PAB.getId()} * {C.getId()}'),
      ('prop_PAC_A', {PAC.getId(): 1, A.getId(): 1}, {PCA.getId(): 1}, f'{kpCA.getId()} * {PAC.getId()} * {A.getId()}'),
      ('prop_PAC_B', {PAC.getId(): 1, B.getId(): 1}, {PCB.getId(): 1}, f'{kpCB.getId()} * {PAC.getId()} * {B.getId()}'),
      ('prop_PAC_C', {PAC.getId(): 1, C.getId(): 1}, {PCC.getId(): 1}, f'{kpCC.getId()} * {PAC.getId()} * {C.getId()}'),
      ('prop_PBA_A', {PBA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {PBA.getId()} * {A.getId()}'),
      ('prop_PBA_B', {PBA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {PBA.getId()} * {B.getId()}'),
      ('prop_PBA_C', {PBA.getId(): 1, C.getId(): 1}, {PAC.getId(): 1}, f'{kpAC.getId()} * {PBA.getId()} * {C.getId()}'),
      ('prop_PBB_A', {PBB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {PBB.getId()} * {A.getId()}'),
      ('prop_PBB_B', {PBB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {PBB.getId()} * {B.getId()}'),
      ('prop_PBB_C', {PBB.getId(): 1, C.getId(): 1}, {PBC.getId(): 1}, f'{kpBC.getId()} * {PBB.getId()} * {C.getId()}'),
      ('prop_PBC_A', {PBC.getId(): 1, A.getId(): 1}, {PCA.getId(): 1}, f'{kpCA.getId()} * {PBC.getId()} * {A.getId()}'),
      ('prop_PBC_B', {PBC.getId(): 1, B.getId(): 1}, {PCB.getId(): 1}, f'{kpCB.getId()} * {PBC.getId()} * {B.getId()}'),
      ('prop_PBC_C', {PBC.getId(): 1, C.getId(): 1}, {PCC.getId(): 1}, f'{kpCC.getId()} * {PBC.getId()} * {C.getId()}'),
      ('prop_PCA_A', {PCA.getId(): 1, A.getId(): 1}, {PAA.getId(): 1}, f'{kpAA.getId()} * {PCA.getId()} * {A.getId()}'),
      ('prop_PCA_B', {PCA.getId(): 1, B.getId(): 1}, {PAB.getId(): 1}, f'{kpAB.getId()} * {PCA.getId()} * {B.getId()}'),
      ('prop_PCA_C', {PCA.getId(): 1, C.getId(): 1}, {PAC.getId(): 1}, f'{kpAC.getId()} * {PCA.getId()} * {C.getId()}'),
      ('prop_PCB_A', {PCB.getId(): 1, A.getId(): 1}, {PBA.getId(): 1}, f'{kpBA.getId()} * {PCB.getId()} * {A.getId()}'),
      ('prop_PCB_B', {PCB.getId(): 1, B.getId(): 1}, {PBB.getId(): 1}, f'{kpBB.getId()} * {PCB.getId()} * {B.getId()}'),
      ('prop_PCB_C', {PCB.getId(): 1, C.getId(): 1}, {PBC.getId(): 1}, f'{kpBC.getId()} * {PCB.getId()} * {C.getId()}'),
      ('prop_PCC_A', {PCC.getId(): 1, A.getId(): 1}, {PCA.getId(): 1}, f'{kpCA.getId()} * {PCC.getId()} * {A.getId()}'),
      ('prop_PCC_B', {PCC.getId(): 1, B.getId(): 1}, {PCB.getId(): 1}, f'{kpCB.getId()} * {PCC.getId()} * {B.getId()}'),
      ('prop_PCC_C', {PCC.getId(): 1, C.getId(): 1}, {PCC.getId(): 1}, f'{kpCC.getId()} * {PCC.getId()} * {C.getId()}'),
  ]

  # Define assignment rules
  create_rule(model, kpBB, formula=f'{kpAA.getId()} / {rAxB.getId()}') # rAxB = kpAA / kpBB
  create_rule(model, kpCC, formula=f'{kpAA.getId()} / {rAxC.getId()}') # rAxC = kpAA / kpCC
  
  create_rule(model, kpAB, formula=f'{kpAA.getId()} / {rAB.getId()}') # rAB = kpAA / kpAB
  create_rule(model, kpAC, formula=f'{kpAA.getId()} / {rAC.getId()}') # rAC = kpAA / kpAC
  
  create_rule(model, kpBA, formula=f'{kpBB.getId()} / {rBA.getId()}') # rBA = kpBB / kpBA
  create_rule(model, kpBC, formula=f'{kpBB.getId()} / {rBC.getId()}') # rBC = kpBB / kpBC
  
  create_rule(model, kpCA, formula=f'{kpCC.getId()} / {rCA.getId()}') # rCA = kpCC / kpCA
  create_rule(model, kpCB, formula=f'{kpCC.getId()} / {rCB.getId()}') # rCB = kpCC / kpCB

  # (reaction_id, reactants_dict, products_dict, kinetic_law)
  print('Creating reactions')
  generated_reactions = []
  for r in reactions:
      reaction = create_reaction(model, r[0], r[1], r[2], r[3])
      generated_reactions.append(reaction)
      
  outputSBML(document, model_filepath)
  
  return model
  