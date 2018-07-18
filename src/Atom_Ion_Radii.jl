# To (somewhat) accurately plot crystal structures, it is (somewhat) useful
# to be able to look-up Atom and Ion approximate radii.
#
# The data contained in this file comes from WebElements.
export element_radius
function element_radius(AE::AbstractString)
  (AN,El,Is)=elementsplit(AE) # Atomic Number, Element symbol, Ionization state
  r=NaN
  if isempty(Is) # no ion specification, return the elemental radius
    if haskey(atomic_radius_calculated,El)
      r=atomic_radius_calculated[El]
    elseif haskey(atomic_radius_empirical,El)
      r=atomic_radius_empirical[El]
    end
  else
    i=isa(parse(Is),Integer)?parse(Is):isa(parse(reverse(Is)),Integer)?parse(reverse(Is)):isa(parse(Is*"1"),Integer)?parse(Is*"1"):0
    y=findfirst(ionic_radius_ionzation_states.==i)
    x=findfirst(ionic_radius_elements.==El)
    x>0&&y>0&&(r=ionic_radius_effective[x,y])
  end
  return isnan(r)?1.:r/100. # default to 1 Å, otherwise convert r from pm to Å
end


# Empirical atomic radii, starting from a data table in the source of
# https://www.webelements.com/periodicity/atomic_radius_empirical/
# extracted via:
# sed -e "s/.*, value:\([^,]\+\),.*, symbol:'\([^']\+\).*/(:\"\2\",\1),/" -e "/null/d"
atomic_radius_empirical=Dict(((:"H",25),(:"Li",145),(:"Be",105),(:"B",85),(:"C",70),(:"N",65),
                             (:"O",60),(:"F",50),(:"Na",180),(:"Mg",150),(:"Al",125),(:"Si",110),
                             (:"P",100),(:"S",100),(:"Cl",100),(:"K",220),(:"Ca",180),(:"Sc",160),
                             (:"Ti",140),(:"V",135),(:"Cr",140),(:"Mn",140),(:"Fe",140),(:"Co",135),
                             (:"Ni",135),(:"Cu",135),(:"Zn",135),(:"Ga",130),(:"Ge",125),(:"As",115),
                             (:"Se",115),(:"Br",115),(:"Rb",235),(:"Sr",200),(:"Y",180),(:"Zr",155),
                             (:"Nb",145),(:"Mo",145),(:"Tc",135),(:"Ru",130),(:"Rh",135),(:"Pd",140),
                             (:"Ag",160),(:"Cd",155),(:"In",155),(:"Sn",145),(:"Sb",145),(:"Te",140),
                             (:"I",140),(:"Cs",260),(:"Ba",215),(:"La",195),(:"Ce",185),(:"Pr",185),
                             (:"Nd",185),(:"Pm",185),(:"Sm",185),(:"Eu",185),(:"Gd",180),(:"Tb",175),
                             (:"Dy",175),(:"Ho",175),(:"Er",175),(:"Tm",175),(:"Yb",175),(:"Lu",175),
                             (:"Hf",155),(:"Ta",145),(:"W",135),(:"Re",135),(:"Os",130),(:"Ir",135),
                             (:"Pt",135),(:"Au",135),(:"Hg",150),(:"Tl",190),(:"Pb",180),(:"Bi",160),
                             (:"Po",190),(:"Ra",215),(:"Ac",195),(:"Th",180),(:"Pa",180),(:"U",175),
                             (:"Np",175),(:"Pu",175),(:"Am",175)));
# Calculated atomic radii, starting from a data table in the source of
# https://www.webelements.com/periodicity/atomic_radius/
atomic_radius_calculated=Dict(((:"H",53),(:"He",31),(:"Li",167),(:"Be",112),(:"B",87),(:"C",67),
                              (:"N",56),(:"O",48),(:"F",42),(:"Ne",38),(:"Na",190),(:"Mg",145),
                              (:"Al",118),(:"Si",111),(:"P",98),(:"S",88),(:"Cl",79),(:"Ar",71),
                              (:"K",243),(:"Ca",194),(:"Sc",184),(:"Ti",176),(:"V",171),(:"Cr",166),
                              (:"Mn",161),(:"Fe",156),(:"Co",152),(:"Ni",149),(:"Cu",145),(:"Zn",142),
                              (:"Ga",136),(:"Ge",125),(:"As",114),(:"Se",103),(:"Br",94),(:"Kr",88),
                              (:"Rb",265),(:"Sr",219),(:"Y",212),(:"Zr",206),(:"Nb",198),(:"Mo",190),
                              (:"Tc",183),(:"Ru",178),(:"Rh",173),(:"Pd",169),(:"Ag",165),(:"Cd",161),
                              (:"In",156),(:"Sn",145),(:"Sb",133),(:"Te",123),(:"I",115),(:"Xe",108),
                              (:"Cs",298),(:"Ba",253),(:"Pr",247),(:"Nd",206),(:"Pm",205),(:"Sm",238),
                              (:"Eu",231),(:"Gd",233),(:"Tb",225),(:"Dy",228),(:"Ho",226),(:"Er",226),
                              (:"Tm",222),(:"Yb",222),(:"Lu",217),(:"Hf",208),(:"Ta",200),(:"W",193),
                              (:"Re",188),(:"Os",185),(:"Ir",180),(:"Pt",177),(:"Au",174),(:"Hg",171),
                              (:"Tl",156),(:"Pb",154),(:"Bi",143),(:"Po",135),(:"At",127),(:"Rn",120)));

# And the effective ionic radii come from
# R. D. Shannon (1976). "Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides".
# Acta Crystallogr A. 32: 751–767. doi:10.1107/S0567739476001551.
# via Wikipedia and can be extracted via
# wget https://en.wikipedia.org/wiki/Ionic_radius -O- | sed -n "/caption><u><i>Effective/,/table/p" | sed -e '/Name/d' -e 's|<sup.*/sup>||' -e '/<a.*\/a>/d' -e 's|(.*)||' -e '/br/d' -e '/table/d' -e '/caption/d' -e 's/[–,−]/-/g' -e 's|<i>[^<]*</i>||g' | sed -r 's|<\/?b>||g'| sed -r 's|<t[h,d]>([\+,-]?[0-9,a-z,A-Z,\.]+)([\+,-]?)\s*<\/t[h,d]>|\2\1|' | sed -e 's|<td></td>|NaN|' | paste -sd " " | sed -r -e 's|<tr>([^<]*)</tr>|\1\n|g' | sed -e "s/^\s*//" -e "s/\s*$//" -e "s/\s/\t/g" > intermediate.txt
#
# to generate the list of elements: cat intermediate.txt | cut -f 2 | sed -r "s/(.*)/\"\1\"/" | paste -sd","
# NOTE: some elements are duplicated due to high-spin/low-spin states. for our current purposes this is fine to ignore
ionic_radius_elements=["H","Li","Be","B","C","N","O","F","Na","Mg","Al","Si",
                       "P","S","Cl","K","Ca","Sc","Ti","V","Cr","Cr","Mn","Mn",
                       "Fe","Fe","Co","Co","Ni","Ni","Cu","Zn","Ga","Ge","As",
                       "Se","Br","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh",
                       "Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba",
                       "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho",
                       "Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",
                       "Au","Hg","Tl","Pb","Bi","Po","At","Fr","Ra","Ac","Th",
                       "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es"];
ionic_radius_ionzation_states=[-3,-2,-1,+1,+2,+3,+4,+5,+6,+7,+8];
ionic_radius_effective=[NaN	NaN	NaN	10	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	76	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	45	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	27	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	16	NaN	NaN	NaN	NaN
                        146	NaN	NaN	NaN	NaN	16	NaN	13	NaN	NaN	NaN
                        NaN	140	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	133	NaN	NaN	NaN	NaN	NaN	NaN	8	NaN
                        NaN	NaN	NaN	102	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	72	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	53.5	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	40	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	44	NaN	38	NaN	NaN	NaN
                        NaN	184	NaN	NaN	NaN	NaN	37	NaN	29	NaN	NaN
                        NaN	NaN	181	NaN	NaN	NaN	NaN	12	NaN	27	NaN
                        NaN	NaN	NaN	138	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	100	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	74.5	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	86	67	60.5	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	79	64	58	54	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	73	61.5	55	49	44	NaN	NaN
                        NaN	NaN	NaN	NaN	80	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	67	58	53	33	25.5	46	NaN
                        NaN	NaN	NaN	NaN	83	64.5	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	61	55	58.5	NaN	25	NaN	NaN
                        NaN	NaN	NaN	NaN	78	64.5	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	65	54.5	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	74.5	61	53	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	69	56	48	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	60	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	77	73	54	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	74	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	62	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	73	NaN	53	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	58	NaN	46	NaN	NaN	NaN
                        NaN	198	NaN	NaN	NaN	NaN	50	NaN	42	NaN	NaN
                        NaN	NaN	196	NaN	NaN	59	NaN	31	NaN	39	NaN
                        NaN	NaN	NaN	152	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	118	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	90	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	72	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	72	68	64	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	69	65	61	59	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	64.5	60	NaN	56	NaN
                        NaN	NaN	NaN	NaN	NaN	68	62	56.5	NaN	38	36
                        NaN	NaN	NaN	NaN	NaN	66.5	60	55	NaN	NaN	NaN
                        NaN	NaN	NaN	59	86	76	61.5	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	115	94	75	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	95	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	80	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	69	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	76	NaN	60	NaN	NaN	NaN
                        NaN	221	NaN	NaN	NaN	NaN	97	NaN	56	NaN	NaN
                        NaN	NaN	220	NaN	NaN	NaN	NaN	95	NaN	53	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	48
                        NaN	NaN	NaN	167	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	135	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	103.2	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	101	87	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	99	85	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	129	98.3	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	97	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	122	95.8	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	117	94.7	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	93.5	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	92.3	76	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	107	91.2	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	90.1	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	89	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	103	88	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	102	86.8	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	86.1	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	71	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	72	68	64	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	66	62	60	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	63	58	55	53	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	63	57.5	54.5	52.5	39
                        NaN	NaN	NaN	NaN	NaN	68	62.5	57	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	80	NaN	62.5	57	NaN	NaN	NaN
                        NaN	NaN	NaN	137	NaN	85	NaN	57	NaN	NaN	NaN
                        NaN	NaN	NaN	119	102	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	150	NaN	88.5	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	119	NaN	77.5	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	103	NaN	76	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	94	NaN	67	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	62	NaN
                        NaN	NaN	NaN	180	NaN	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	148	NaN	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	112	NaN	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	NaN	94	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	104	90	78	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	102.5	89	76	73	NaN	NaN
                        NaN	NaN	NaN	NaN	110	101	87	75	72	71	NaN
                        NaN	NaN	NaN	NaN	NaN	100	86	74	71	NaN	NaN
                        NaN	NaN	NaN	NaN	126	97.5	85	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	97	85	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	96	83	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	95	82.1	NaN	NaN	NaN	NaN
                        NaN	NaN	NaN	NaN	NaN	83.5	NaN	NaN	NaN	NaN	NaN];
