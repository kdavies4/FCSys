(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[         0,          0]
NotebookDataLength[     62671,       1912]
NotebookOptionsPosition[     60625,       1848]
NotebookOutlinePosition[     60962,       1863]
CellTagsIndexPosition[     60919,       1860]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell[BoxData[
 RowBox[{"(*", " ",
  RowBox[{
   RowBox[{
   "Solve", " ", "the", " ", "systems", " ", "of", " ", "units", " ",
    RowBox[{"(",
     RowBox[{
     "m", ",", " ", "kg", ",", " ", "s", ",", " ", "A", ",", " ", "K", ",",
      " ",
      RowBox[{"and", " ", "mol"}]}], ")"}], " ", "for", " ", "the", " ",
    "Rydberg", " ", "constant", " ",
    RowBox[{"(", "Rinf", ")"}]}], ",", " ",
   RowBox[{"speed", " ", "of", " ", "light", " ", "in", " ", "vacuum", " ",
    RowBox[{"(", "c0", ")"}]}], ",", " ",
   RowBox[{"von", " ", "Klitzing", " ", "constant", " ",
    RowBox[{"(", "kvK", ")"}]}], ",", " ",
   RowBox[{"and", " ", "Josephson", " ", "constant", " ",
    RowBox[{"(", "kJ", ")"}], " ", "by", " ", "normalizing", " ",
    RowBox[{"(",
     RowBox[{"setting", " ", "to", " ", "unity"}], ")"}], " ", "four", " ",
    "of", " ", "the", " ",
    RowBox[{"units", ".", "  ", "The"}], " ", "Faraday", " ", "and", " ",
    "gas", " ", "constants", " ", "are", " ", "always", " ",
    RowBox[{"normalized", ".", " ", "\[IndentingNewLine]", "Kevin"}], " ",
    RowBox[{"L", ".", " ", "Davies"}], "\[IndentingNewLine]",
    RowBox[{
     RowBox[{"6", "/", "4"}], "/", "2012"}]}]}], " ", "*)"}]], "Input",
 CellChangeTimes->{{3.547829105305603*^9, 3.547829199912115*^9}, {
   3.547829238291985*^9, 3.547829446574378*^9}, {3.547829671002751*^9,
   3.547829691974365*^9}, {3.547829737004084*^9, 3.547829750522901*^9}, {
   3.547829782045839*^9, 3.547829782283955*^9}, {3.547830336932003*^9,
   3.547830337464282*^9}, {3.547853917121263*^9, 3.547853936550551*^9}, {
   3.547853970138355*^9, 3.547854286926183*^9}, 3.547856973863179*^9, {
   3.548893586055835*^9, 3.548893636446683*^9}, {3.548893668869068*^9,
   3.548893673554918*^9}}],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"x", "=",
  RowBox[{"{",
   RowBox[{"10973731.568539", ",", "299792458", ",",
    RowBox[{"483597.870*^9", "*",
     RowBox[{
      RowBox[{"(",
       RowBox[{"Pi", "*", "1*^-7"}], ")"}], "^", "2"}]}], ",", "96485.3365",
    ",", "8.3144621", ",", "25812.8074434"}],
   "}"}]}], "\[IndentingNewLine]",
 RowBox[{"a", "=",
  RowBox[{"{",
   RowBox[{
    RowBox[{
     RowBox[{"x", "[",
      RowBox[{"[", "1", "]"}], "]"}], "*",
     RowBox[{"x", "[",
      RowBox[{"[", "4", "]"}], "]"}], "*",
     RowBox[{
      RowBox[{
       RowBox[{"x", "[",
        RowBox[{"[", "6", "]"}], "]"}], "^", "2"}], "/",
      RowBox[{"(",
       RowBox[{
        RowBox[{"x", "[",
         RowBox[{"[", "3", "]"}], "]"}], "*",
        RowBox[{
         RowBox[{"x", "[",
          RowBox[{"[", "2", "]"}], "]"}], "^", "3"}]}], ")"}]}]}], ",",
    RowBox[{"1", "/",
     RowBox[{"x", "[",
      RowBox[{"[", "4", "]"}], "]"}]}], ",",
    RowBox[{
     RowBox[{"x", "[",
      RowBox[{"[", "3", "]"}], "]"}], "*",
     RowBox[{"x", "[",
      RowBox[{"[", "5", "]"}], "]"}], "*",
     RowBox[{
      RowBox[{
       RowBox[{"x", "[",
        RowBox[{"[", "2", "]"}], "]"}], "^", "3"}], "/",
      RowBox[{"(",
       RowBox[{
        RowBox[{"x", "[",
         RowBox[{"[", "1", "]"}], "]"}], "*",
        RowBox[{"x", "[",
         RowBox[{"[", "4", "]"}], "]"}], "*",
        RowBox[{
         RowBox[{"x", "[",
          RowBox[{"[", "6", "]"}], "]"}], "^", "2"}]}], ")"}]}]}]}],
   "}"}]}], "\[IndentingNewLine]",
 RowBox[{"system1", ":=",
  RowBox[{
   RowBox[{
    RowBox[{"A", "*",
     RowBox[{"s", "/", "mol"}]}], "\[Equal]",
    RowBox[{"a", "[",
     RowBox[{"[", "2", "]"}], "]"}]}], "&&",
   RowBox[{
    RowBox[{"kelvin", "*",
     RowBox[{
      RowBox[{"s", "^", "2"}], "/",
      RowBox[{"m", "^", "2"}]}]}], "\[Equal]",
    RowBox[{"a", "[",
     RowBox[{"[", "3", "]"}], "]"}]}]}]}]}], "Input",
 CellChangeTimes->{{3.547856869845613*^9, 3.547856873097409*^9}, {
   3.547856925074122*^9, 3.547856925593602*^9}, {3.547861148466646*^9,
   3.547861171789515*^9}, {3.547861218636286*^9, 3.547861219092502*^9}, {
   3.547861314356918*^9, 3.547861376623362*^9}, {3.54786163899545*^9,
   3.547861657784416*^9}, {3.547864080114068*^9, 3.54786411010159*^9}, {
   3.547864157562982*^9, 3.547864181398806*^9}, 3.548893906119694*^9, {
   3.54894395080405*^9, 3.548944018246586*^9}, {3.548960587828049*^9,
   3.548960611877159*^9}, {3.55518428398876*^9, 3.555184293166202*^9},
   3.555186691445013*^9, {3.555186836668636*^9, 3.555186856899824*^9}, {
   3.555187366887889*^9, 3.555187367058878*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
  "1.0973731568539`*^7", ",", "299792458", ",", "47.72919666109439`", ",",
   "96485.3365`", ",", "8.3144621`", ",", "25812.8074434`"}], "}"}]], "Output",\

 CellChangeTimes->{{3.547856874886548*^9, 3.547856926351463*^9},
   3.547857240084261*^9, 3.547861671727941*^9, 3.5478642764493*^9,
   3.54786658652126*^9, 3.548893914306848*^9, {3.548944021174859*^9,
   3.548944038401247*^9}, 3.548960403628301*^9, 3.55518383872703*^9,
   3.555184306835867*^9, 3.555184364924438*^9, 3.555186743750541*^9, {
   3.555186824971502*^9, 3.555186841225712*^9}, 3.555186878678313*^9,
   3.555187254430433*^9, 3.555187424959131*^9, 3.555187467970709*^9}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
  "5.485799092185156`*^-7", ",", "0.000010364269186126536`", ",",
   "1.5156337226867167`*^7"}], "}"}]], "Output",
 CellChangeTimes->{{3.547856874886548*^9, 3.547856926351463*^9},
   3.547857240084261*^9, 3.547861671727941*^9, 3.5478642764493*^9,
   3.54786658652126*^9, 3.548893914306848*^9, {3.548944021174859*^9,
   3.548944038401247*^9}, 3.548960403628301*^9, 3.55518383872703*^9,
   3.555184306835867*^9, 3.555184364924438*^9, 3.555186743750541*^9, {
   3.555186824971502*^9, 3.555186841225712*^9}, 3.555186878678313*^9,
   3.555187254430433*^9, 3.555187424959131*^9, 3.555187467974754*^9}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ",
   RowBox[{
    RowBox[{
     RowBox[{
      RowBox[{"**", "8"}], "/", "28"}], "/", "12"}], ":", " ",
    RowBox[{
    "These", " ", "values", " ", "do", " ", "not", " ", "seem", " ", "to",
     " ", "match", " ", "those", " ", "of", " ", "Modelica", " ",
     RowBox[{
      RowBox[{"(",
       RowBox[{"FCSys", ".", "Units"}], ")"}], "."}]}]}], " ", "*)"}],
  "\[IndentingNewLine]",
  RowBox[{
   RowBox[{
    RowBox[{"Solve", "[",
     RowBox[{"system1", ",",
      RowBox[{"{",
       RowBox[{"mol", ",", "s"}], "}"}]}], "]"}], "/.",
    RowBox[{"{",
     RowBox[{
      RowBox[{"A", "\[Rule]", "1"}], ",",
      RowBox[{"kelvin", "\[Rule]", "1"}], ",",
      RowBox[{"m", "\[Rule]", "1"}], ",",
      RowBox[{"kg", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
   RowBox[{
    RowBox[{"Solve", "[",
     RowBox[{"system1", ",",
      RowBox[{"{",
       RowBox[{"m", ",", "s"}], "}"}]}], "]"}], "/.",
    RowBox[{"{",
     RowBox[{
      RowBox[{"A", "\[Rule]", "1"}], ",",
      RowBox[{"kelvin", "\[Rule]", "1"}], ",",
      RowBox[{"mol", "\[Rule]", "1"}], ",",
      RowBox[{"kg", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
   RowBox[{
    RowBox[{"Solve", "[",
     RowBox[{"system1", ",",
      RowBox[{"{",
       RowBox[{"m", ",", "mol"}], "}"}]}], "]"}], "/.",
    RowBox[{"{",
     RowBox[{
      RowBox[{"A", "\[Rule]", "1"}], ",",
      RowBox[{"kelvin", "\[Rule]", "1"}], ",",
      RowBox[{"s", "\[Rule]", "1"}], ",",
      RowBox[{"kg", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
   RowBox[{
    RowBox[{"Solve", "[",
     RowBox[{"system1", ",",
      RowBox[{"{",
       RowBox[{"kelvin", ",", "s"}], "}"}]}], "]"}], "/.",
    RowBox[{"{",
     RowBox[{
      RowBox[{"A", "\[Rule]", "1"}], ",",
      RowBox[{"m", "\[Rule]", "1"}], ",",
      RowBox[{"mol", "\[Rule]", "1"}], ",",
      RowBox[{"kg", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
   RowBox[{
    RowBox[{"Solve", "[",
     RowBox[{"system1", ",",
      RowBox[{"{",
       RowBox[{"kelvin", ",", "mol"}], "}"}]}], "]"}], "/.",
    RowBox[{"{",
     RowBox[{
      RowBox[{"A", "\[Rule]", "1"}], ",",
      RowBox[{"m", "\[Rule]", "1"}], ",",
      RowBox[{"s", "\[Rule]", "1"}], ",",
      RowBox[{"kg", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
   RowBox[{
    RowBox[{"Solve", "[",
     RowBox[{"system1", ",",
      RowBox[{"{",
       RowBox[{"A", ",", "s"}], "}"}]}], "]"}], "/.",
    RowBox[{"{",
     RowBox[{
      RowBox[{"kelvin", "\[Rule]", "1"}], ",",
      RowBox[{"m", "\[Rule]", "1"}], ",",
      RowBox[{"mol", "\[Rule]", "1"}], ",",
      RowBox[{"kg", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
   RowBox[{
    RowBox[{"Solve", "[",
     RowBox[{"system1", ",",
      RowBox[{"{",
       RowBox[{"m", ",", "A"}], "}"}]}], "]"}], "/.",
    RowBox[{"{",
     RowBox[{
      RowBox[{"kelvin", "\[Rule]", "1"}], ",",
      RowBox[{"mol", "\[Rule]", "1"}], ",",
      RowBox[{"s", "\[Rule]", "1"}], ",",
      RowBox[{"kg", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
   RowBox[{
    RowBox[{"Solve", "[",
     RowBox[{"system1", ",",
      RowBox[{"{",
       RowBox[{"kelvin", ",", "A"}], "}"}]}], "]"}], "/.",
    RowBox[{"{",
     RowBox[{
      RowBox[{"m", "\[Rule]", "1"}], ",",
      RowBox[{"mol", "\[Rule]", "1"}], ",",
      RowBox[{"s", "\[Rule]", "1"}], ",",
      RowBox[{"kg", "\[Rule]", "1"}]}], "}"}]}]}]}]], "Input",
 CellChangeTimes->{{3.547856899047643*^9, 3.547856909055331*^9}, {
  3.547857129850629*^9, 3.547857143221439*^9}, {3.548893714099826*^9,
  3.54889381261119*^9}, {3.548893973528711*^9, 3.548894080953529*^9}, {
  3.555185862406578*^9, 3.555185903016705*^9}, {3.555186675698845*^9,
  3.555186702132998*^9}, {3.555186734525444*^9, 3.555186735036974*^9}, {
  3.55518683674189*^9, 3.555186857028001*^9}, {3.555187335003345*^9,
  3.555187358546286*^9}, {3.555187443365893*^9, 3.555187449863292*^9}, {
  3.555187486362729*^9, 3.55518749092181*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Solve", "::", "ratnz"}], "MessageName"],
  RowBox[{
  ":", " "}], "\<\"Solve was unable to solve the system with inexact \
coefficients. The answer was obtained by solving a corresponding exact system \
and numericizing the result. \\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", \
ButtonStyle->\\\"Link\\\", ButtonFrame->None, \
ButtonData:>\\\"paclet:ref/Solve\\\", ButtonNote -> \
\\\"Solve::ratnz\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.555187468110226*^9}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"mol", "\[Rule]",
      RowBox[{"-", "3.7562842188025904`*^8"}]}], ",",
     RowBox[{"s", "\[Rule]",
      RowBox[{"-", "3893.1140783269075`"}]}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"mol", "\[Rule]", "3.7562842188025904`*^8"}], ",",
     RowBox[{"s", "\[Rule]", "3893.1140783269075`"}]}], "}"}]}],
  "}"}]], "Output",
 CellChangeTimes->{
  3.555186743857704*^9, {3.555186825240059*^9, 3.555186841331698*^9},
   3.555186878954844*^9, 3.555187254590753*^9, {3.5551874250718*^9,
   3.555187468113133*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Solve", "::", "ratnz"}], "MessageName"],
  RowBox[{
  ":", " "}], "\<\"Solve was unable to solve the system with inexact \
coefficients. The answer was obtained by solving a corresponding exact system \
and numericizing the result. \\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", \
ButtonStyle->\\\"Link\\\", ButtonFrame->None, \
ButtonData:>\\\"paclet:ref/Solve\\\", ButtonNote -> \
\\\"Solve::ratnz\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.555187468234463*^9}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"m", "\[Rule]",
      RowBox[{"-", "2.6622053650636027`*^-9"}]}], ",",
     RowBox[{"s", "\[Rule]", "0.000010364269186126537`"}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"m", "\[Rule]", "2.6622053650636027`*^-9"}], ",",
     RowBox[{"s", "\[Rule]", "0.000010364269186126537`"}]}], "}"}]}],
  "}"}]], "Output",
 CellChangeTimes->{
  3.555186743857704*^9, {3.555186825240059*^9, 3.555186841331698*^9},
   3.555186878954844*^9, 3.555187254590753*^9, {3.5551874250718*^9,
   3.555187468237498*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Solve", "::", "ratnz"}], "MessageName"],
  RowBox[{
  ":", " "}], "\<\"Solve was unable to solve the system with inexact \
coefficients. The answer was obtained by solving a corresponding exact system \
and numericizing the result. \\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", \
ButtonStyle->\\\"Link\\\", ButtonFrame->None, \
ButtonData:>\\\"paclet:ref/Solve\\\", ButtonNote -> \
\\\"Solve::ratnz\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.555187468239252*^9}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"m", "\[Rule]",
      RowBox[{"-", "0.0002568637804802671`"}]}], ",",
     RowBox[{"mol", "\[Rule]", "96485.3365`"}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"m", "\[Rule]", "0.0002568637804802671`"}], ",",
     RowBox[{"mol", "\[Rule]", "96485.3365`"}]}], "}"}]}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186743857704*^9, {3.555186825240059*^9, 3.555186841331698*^9},
   3.555186878954844*^9, 3.555187254590753*^9, {3.5551874250718*^9,
   3.555187468241012*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Solve", "::", "ratnz"}], "MessageName"],
  RowBox[{
  ":", " "}], "\<\"Solve was unable to solve the system with inexact \
coefficients. The answer was obtained by solving a corresponding exact system \
and numericizing the result. \\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", \
ButtonStyle->\\\"Link\\\", ButtonFrame->None, \
ButtonData:>\\\"paclet:ref/Solve\\\", ButtonNote -> \
\\\"Solve::ratnz\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.555187468242655*^9}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{"{",
   RowBox[{
    RowBox[{"kelvin", "\[Rule]", "1.4109671132425387`*^17"}], ",",
    RowBox[{"s", "\[Rule]", "0.000010364269186126537`"}]}], "}"}],
  "}"}]], "Output",
 CellChangeTimes->{
  3.555186743857704*^9, {3.555186825240059*^9, 3.555186841331698*^9},
   3.555186878954844*^9, 3.555187254590753*^9, {3.5551874250718*^9,
   3.555187468244318*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Solve", "::", "ratnz"}], "MessageName"],
  RowBox[{
  ":", " "}], "\<\"Solve was unable to solve the system with inexact \
coefficients. The answer was obtained by solving a corresponding exact system \
and numericizing the result. \\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", \
ButtonStyle->\\\"Link\\\", ButtonFrame->None, \
ButtonData:>\\\"paclet:ref/Solve\\\", ButtonNote -> \
\\\"Solve::ratnz\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.555187468245874*^9}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{"{",
   RowBox[{
    RowBox[{"kelvin", "\[Rule]", "1.515633722686717`*^7"}], ",",
    RowBox[{"mol", "\[Rule]", "96485.3365`"}]}], "}"}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186743857704*^9, {3.555186825240059*^9, 3.555186841331698*^9},
   3.555186878954844*^9, 3.555187254590753*^9, {3.5551874250718*^9,
   3.555187468247547*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Solve", "::", "ratnz"}], "MessageName"],
  RowBox[{
  ":", " "}], "\<\"Solve was unable to solve the system with inexact \
coefficients. The answer was obtained by solving a corresponding exact system \
and numericizing the result. \\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", \
ButtonStyle->\\\"Link\\\", ButtonFrame->None, \
ButtonData:>\\\"paclet:ref/Solve\\\", ButtonNote -> \
\\\"Solve::ratnz\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.55518746824906*^9}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"A", "\[Rule]",
      RowBox[{"-", "2.6622053650636027`*^-9"}]}], ",",
     RowBox[{"s", "\[Rule]",
      RowBox[{"-", "3893.1140783269075`"}]}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"A", "\[Rule]", "2.6622053650636027`*^-9"}], ",",
     RowBox[{"s", "\[Rule]", "3893.1140783269075`"}]}], "}"}]}],
  "}"}]], "Output",
 CellChangeTimes->{
  3.555186743857704*^9, {3.555186825240059*^9, 3.555186841331698*^9},
   3.555186878954844*^9, 3.555187254590753*^9, {3.5551874250718*^9,
   3.555187468250829*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Solve", "::", "ratnz"}], "MessageName"],
  RowBox[{
  ":", " "}], "\<\"Solve was unable to solve the system with inexact \
coefficients. The answer was obtained by solving a corresponding exact system \
and numericizing the result. \\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", \
ButtonStyle->\\\"Link\\\", ButtonFrame->None, \
ButtonData:>\\\"paclet:ref/Solve\\\", ButtonNote -> \
\\\"Solve::ratnz\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.555187468252565*^9}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"m", "\[Rule]",
      RowBox[{"-", "0.0002568637804802671`"}]}], ",",
     RowBox[{"A", "\[Rule]", "0.000010364269186126537`"}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"m", "\[Rule]", "0.0002568637804802671`"}], ",",
     RowBox[{"A", "\[Rule]", "0.000010364269186126537`"}]}], "}"}]}],
  "}"}]], "Output",
 CellChangeTimes->{
  3.555186743857704*^9, {3.555186825240059*^9, 3.555186841331698*^9},
   3.555186878954844*^9, 3.555187254590753*^9, {3.5551874250718*^9,
   3.555187468254332*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{"{",
   RowBox[{
    RowBox[{"kelvin", "\[Rule]", "1.5156337226867167`*^7"}], ",",
    RowBox[{"A", "\[Rule]", "0.000010364269186126536`"}]}], "}"}],
  "}"}]], "Output",
 CellChangeTimes->{
  3.555186743857704*^9, {3.555186825240059*^9, 3.555186841331698*^9},
   3.555186878954844*^9, 3.555187254590753*^9, {3.5551874250718*^9,
   3.555187468255966*^9}}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"x", ":=",
  RowBox[{"{",
   RowBox[{"x1", ",", "x2", ",", "x3", ",", "x4", ",", "x5", ",", "x6"}],
   "}"}]}], "\[IndentingNewLine]",
 RowBox[{"a", "=",
  RowBox[{"{",
   RowBox[{
    RowBox[{
     RowBox[{"x", "[",
      RowBox[{"[", "1", "]"}], "]"}], "*",
     RowBox[{"x", "[",
      RowBox[{"[", "4", "]"}], "]"}], "*",
     RowBox[{
      RowBox[{
       RowBox[{"x", "[",
        RowBox[{"[", "6", "]"}], "]"}], "^", "2"}], "/",
      RowBox[{"(",
       RowBox[{
        RowBox[{"x", "[",
         RowBox[{"[", "3", "]"}], "]"}], "*",
        RowBox[{
         RowBox[{"x", "[",
          RowBox[{"[", "2", "]"}], "]"}], "^", "3"}]}], ")"}]}]}], ",",
    RowBox[{"1", "/",
     RowBox[{"x", "[",
      RowBox[{"[", "4", "]"}], "]"}]}], ",",
    RowBox[{
     RowBox[{"x", "[",
      RowBox[{"[", "3", "]"}], "]"}], "*",
     RowBox[{"x", "[",
      RowBox[{"[", "5", "]"}], "]"}], "*",
     RowBox[{
      RowBox[{
       RowBox[{"x", "[",
        RowBox[{"[", "2", "]"}], "]"}], "^", "3"}], "/",
      RowBox[{"(",
       RowBox[{
        RowBox[{"x", "[",
         RowBox[{"[", "1", "]"}], "]"}], "*",
        RowBox[{"x", "[",
         RowBox[{"[", "4", "]"}], "]"}], "*",
        RowBox[{
         RowBox[{"x", "[",
          RowBox[{"[", "6", "]"}], "]"}], "^", "2"}]}], ")"}]}]}]}],
   "}"}]}], "\[IndentingNewLine]",
 RowBox[{"b", "=",
  RowBox[{"{",
   RowBox[{
    RowBox[{"x", "[",
     RowBox[{"[", "1", "]"}], "]"}], ",",
    RowBox[{"x", "[",
     RowBox[{"[", "2", "]"}], "]"}], ",",
    RowBox[{
     RowBox[{
      RowBox[{"(",
       RowBox[{
        RowBox[{"x", "[",
         RowBox[{"[", "3", "]"}], "]"}], "*",
        RowBox[{
         RowBox[{"x", "[",
          RowBox[{"[", "2", "]"}], "]"}], "^", "3"}]}], ")"}], "^", "2"}],
     "/",
     RowBox[{"(",
      RowBox[{
       RowBox[{
        RowBox[{"x", "[",
         RowBox[{"[", "1", "]"}], "]"}], "^", "2"}], "*",
       RowBox[{
        RowBox[{"x", "[",
         RowBox[{"[", "6", "]"}], "]"}], "^", "3"}]}], ")"}]}]}],
   "}"}]}], "\[IndentingNewLine]",
 RowBox[{"system1", ":=",
  RowBox[{
   RowBox[{
    RowBox[{"A", "*",
     RowBox[{"s", "/", "mol"}]}], "\[Equal]",
    RowBox[{"a", "[",
     RowBox[{"[", "2", "]"}], "]"}]}], "&&",
   RowBox[{
    RowBox[{"kelvin", "*",
     RowBox[{
      RowBox[{"s", "^", "2"}], "/",
      RowBox[{"m", "^", "2"}]}]}], "\[Equal]",
    RowBox[{"a", "[",
     RowBox[{"[", "3", "]"}], "]"}]}]}]}], "\n",
 RowBox[{"system2", ":=",
  RowBox[{
   RowBox[{
    RowBox[{"mol", "/", "kg"}], "\[Equal]",
    RowBox[{"a", "[",
     RowBox[{"[", "1", "]"}], "]"}]}], "&&",
   RowBox[{
    RowBox[{"A", "*",
     RowBox[{"s", "/", "mol"}]}], "==",
    RowBox[{"a", "[",
     RowBox[{"[", "2", "]"}], "]"}]}], "&&",
   RowBox[{
    RowBox[{"kelvin", "*",
     RowBox[{
      RowBox[{"s", "^", "2"}], "/",
      RowBox[{"m", "^", "2"}]}]}], "\[Equal]",
    RowBox[{"a", "[",
     RowBox[{"[", "3", "]"}], "]"}]}], "&&",
   RowBox[{"Rinf", "==",
    RowBox[{
     RowBox[{"b", "[",
      RowBox[{"[", "1", "]"}], "]"}], "/", "m"}]}], "&&",
   RowBox[{"c0", "==",
    RowBox[{
     RowBox[{"b", "[",
      RowBox[{"[", "2", "]"}], "]"}], "*",
     RowBox[{"m", "/", "s"}]}]}], "&&",
   RowBox[{"kvK", "==",
    RowBox[{
     RowBox[{"b", "[",
      RowBox[{"[", "3", "]"}], "]"}], "*",
     RowBox[{
      RowBox[{"m", "^", "2"}], "/",
      RowBox[{"(",
       RowBox[{"kg", "*", "s"}], ")"}]}]}]}]}]}]}], "Input",
 CellChangeTimes->{{3.547850495811396*^9, 3.547850542044922*^9}, {
   3.547850575680254*^9, 3.547850643340001*^9}, 3.547851284685791*^9, {
   3.547852419389999*^9, 3.547852428544582*^9}, {3.547852495773295*^9,
   3.547852499947832*^9}, {3.547855568138406*^9, 3.547855570926971*^9}, {
   3.547855630835598*^9, 3.547855631031077*^9}, {3.54785568965465*^9,
   3.547855692896102*^9}, {3.547855749073499*^9, 3.547855809975341*^9}, {
   3.547855846765299*^9, 3.547855924887463*^9}, {3.547855958920451*^9,
   3.547855966075408*^9}, {3.547856072243744*^9, 3.547856447349873*^9}, {
   3.547856530332584*^9, 3.547856579839943*^9}, {3.547856674912238*^9,
   3.547856699409577*^9}, {3.547856754745164*^9, 3.5478567966689*^9},
   3.547856943615778*^9, 3.548893946059728*^9, {3.555186706391722*^9,
   3.555186737714378*^9}, {3.555186837145045*^9, 3.555186855331279*^9},
   3.555187176551426*^9, {3.555187324579416*^9, 3.555187325674816*^9}, {
   3.555187401970386*^9, 3.555187402159678*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   FractionBox[
    RowBox[{"x1", " ", "x4", " ",
     SuperscriptBox["x6", "2"]}],
    RowBox[{
     SuperscriptBox["x2", "3"], " ", "x3"}]], ",",
   FractionBox["1", "x4"], ",",
   FractionBox[
    RowBox[{
     SuperscriptBox["x2", "3"], " ", "x3", " ", "x5"}],
    RowBox[{"x1", " ", "x4", " ",
     SuperscriptBox["x6", "2"]}]]}], "}"}]], "Output",
 CellChangeTimes->{
  3.54785051103632*^9, 3.547850544841611*^9, {3.547850608376215*^9,
   3.547850644010213*^9}, 3.547851253198036*^9, 3.547851285329707*^9,
   3.547851374380715*^9, 3.547852315132664*^9, 3.547852433191671*^9,
   3.547852502698827*^9, {3.547856574621933*^9, 3.547856580650277*^9}, {
   3.547856675582294*^9, 3.54785670080875*^9}, {3.547856772838318*^9,
   3.54785679793409*^9}, {3.547856938163276*^9, 3.547856944482344*^9},
   3.547866453710085*^9, 3.548944039147827*^9, 3.555183839532253*^9,
   3.555184307024344*^9, 3.555184365285933*^9, 3.555186744179541*^9, {
   3.555186825431161*^9, 3.555186841465905*^9}, 3.55518687901754*^9,
   3.555187254629744*^9, 3.555187425206324*^9, {3.555187459780616*^9,
   3.555187468271131*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{"x1", ",", "x2", ",",
   FractionBox[
    RowBox[{
     SuperscriptBox["x2", "6"], " ",
     SuperscriptBox["x3", "2"]}],
    RowBox[{
     SuperscriptBox["x1", "2"], " ",
     SuperscriptBox["x6", "3"]}]]}], "}"}]], "Output",
 CellChangeTimes->{
  3.54785051103632*^9, 3.547850544841611*^9, {3.547850608376215*^9,
   3.547850644010213*^9}, 3.547851253198036*^9, 3.547851285329707*^9,
   3.547851374380715*^9, 3.547852315132664*^9, 3.547852433191671*^9,
   3.547852502698827*^9, {3.547856574621933*^9, 3.547856580650277*^9}, {
   3.547856675582294*^9, 3.54785670080875*^9}, {3.547856772838318*^9,
   3.54785679793409*^9}, {3.547856938163276*^9, 3.547856944482344*^9},
   3.547866453710085*^9, 3.548944039147827*^9, 3.555183839532253*^9,
   3.555184307024344*^9, 3.555184365285933*^9, 3.555186744179541*^9, {
   3.555186825431161*^9, 3.555186841465905*^9}, 3.55518687901754*^9,
   3.555187254629744*^9, 3.555187425206324*^9, {3.555187459780616*^9,
   3.555187468273237*^9}}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{"Rinf", ",", "c0", ",", "kvK", ",", "m", ",", "kg", ",", "mol"}],
      "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"kelvin", "\[Rule]", "1"}], ",",
    RowBox[{"A", "\[Rule]", "1"}], ",",
    RowBox[{"s", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{"Rinf", ",", "c0", ",", "kvK", ",", "m", ",", "mol", ",", "s"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"kelvin", "\[Rule]", "1"}], ",",
    RowBox[{"kg", "\[Rule]", "1"}], ",",
    RowBox[{"A", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{"Rinf", ",", "c0", ",", "kvK", ",", "m", ",", "A", ",", "mol"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"kelvin", "\[Rule]", "1"}], ",",
    RowBox[{"kg", "\[Rule]", "1"}], ",",
    RowBox[{"s", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{"Rinf", ",", "c0", ",", "kvK", ",", "m", ",", "kg", ",", "s"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"kelvin", "\[Rule]", "1"}], ",",
    RowBox[{"mol", "\[Rule]", "1"}], ",",
    RowBox[{"A", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{"Rinf", ",", "c0", ",", "kvK", ",", "m", ",", "A", ",", "kg"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"kelvin", "\[Rule]", "1"}], ",",
    RowBox[{"mol", "\[Rule]", "1"}], ",",
    RowBox[{"s", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{"Rinf", ",", "c0", ",", "kvK", ",", "kg", ",", "mol", ",", "s"}],
      "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"m", "\[Rule]", "1"}], ",",
    RowBox[{"A", "\[Rule]", "1"}], ",",
    RowBox[{"kelvin", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{
     "Rinf", ",", "c0", ",", "kvK", ",", "kelvin", ",", "kg", ",", "mol"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"m", "\[Rule]", "1"}], ",",
    RowBox[{"A", "\[Rule]", "1"}], ",",
    RowBox[{"s", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{
     "Rinf", ",", "c0", ",", "kvK", ",", "kelvin", ",", "kg", ",", "s"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"m", "\[Rule]", "1"}], ",",
    RowBox[{"mol", "\[Rule]", "1"}], ",",
    RowBox[{"A", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{"Rinf", ",", "c0", ",", "kvK", ",", "A", ",", "kg", ",", "s"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"m", "\[Rule]", "1"}], ",",
    RowBox[{"mol", "\[Rule]", "1"}], ",",
    RowBox[{"kelvin", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{
     "Rinf", ",", "c0", ",", "kvK", ",", "kelvin", ",", "mol", ",", "s"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"m", "\[Rule]", "1"}], ",",
    RowBox[{"kg", "\[Rule]", "1"}], ",",
    RowBox[{"A", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{"Rinf", ",", "c0", ",", "kvK", ",", "A", ",", "mol", ",", "s"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"m", "\[Rule]", "1"}], ",",
    RowBox[{"kg", "\[Rule]", "1"}], ",",
    RowBox[{"kelvin", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{
     "Rinf", ",", "c0", ",", "kvK", ",", "kelvin", ",", "A", ",", "mol"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"m", "\[Rule]", "1"}], ",",
    RowBox[{"kg", "\[Rule]", "1"}], ",",
    RowBox[{"s", "\[Rule]", "1"}]}], "}"}]}], "\[IndentingNewLine]",
 RowBox[{
  RowBox[{"Solve", "[",
   RowBox[{"system2", ",",
    RowBox[{"{",
     RowBox[{
     "Rinf", ",", "c0", ",", "kvK", ",", "kelvin", ",", "A", ",", "kg"}],
     "}"}]}], "]"}], "/.",
  RowBox[{"{",
   RowBox[{
    RowBox[{"m", "\[Rule]", "1"}], ",",
    RowBox[{"mol", "\[Rule]", "1"}], ",",
    RowBox[{"s", "\[Rule]", "1"}]}], "}"}]}]}], "Input",
 CellChangeTimes->{{3.547856603268654*^9, 3.547856644852943*^9},
   3.54785722642976*^9, {3.555186719274813*^9, 3.555186728092909*^9}, {
   3.555186837206422*^9, 3.555186856346681*^9}, {3.555187074796356*^9,
   3.555187105707022*^9}, {3.555187189169308*^9, 3.555187202621898*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}],
        RowBox[{
         SqrtBox["x4"], " ", "x6"}]]}]}], ",",
     RowBox[{"c0", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SqrtBox["x2"], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{"x4", " ", "x6"}], "x5"]}], ",",
     RowBox[{"m", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kg", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}],
       RowBox[{"x1", " ",
        SuperscriptBox["x6", "2"]}]]}], ",",
     RowBox[{"mol", "\[Rule]", "x4"}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}],
       RowBox[{
        SqrtBox["x4"], " ", "x6"}]]}], ",",
     RowBox[{"c0", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SqrtBox["x2"], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{"x4", " ", "x6"}], "x5"]}], ",",
     RowBox[{"m", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kg", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}],
       RowBox[{"x1", " ",
        SuperscriptBox["x6", "2"]}]]}], ",",
     RowBox[{"mol", "\[Rule]", "x4"}]}], "}"}]}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468440003*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"9", "/", "2"}]], " ",
         SuperscriptBox["x3",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x5"]}],
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ",
         SuperscriptBox["x6", "3"]}]]}]}], ",",
     RowBox[{"c0", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SqrtBox["x2"], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{"x4", " ", "x6"}], "x5"]}], ",",
     RowBox[{"m", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x1",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x4"], " ",
         SuperscriptBox["x6", "3"]}],
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"9", "/", "2"}]], " ",
         SuperscriptBox["x3",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"mol", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}], ",",
     RowBox[{"s", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"9", "/", "2"}]], " ",
        SuperscriptBox["x3",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x5"]}],
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ",
        SuperscriptBox["x6", "3"]}]]}], ",",
     RowBox[{"c0", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SqrtBox["x2"], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{"x4", " ", "x6"}], "x5"]}], ",",
     RowBox[{"m", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x1",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x4"], " ",
        SuperscriptBox["x6", "3"]}],
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"9", "/", "2"}]], " ",
        SuperscriptBox["x3",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"mol", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}], ",",
     RowBox[{"s", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}]}], "}"}]}], "}"}]], "Output",\

 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468445332*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}],
        RowBox[{
         SqrtBox["x4"], " ", "x6"}]]}]}], ",",
     RowBox[{"c0", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SqrtBox["x2"], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3", " ", "x4"}],
       RowBox[{"x1", " ", "x5", " ", "x6"}]]}], ",",
     RowBox[{"m", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"A", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}], ",",
     RowBox[{"mol", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}],
       RowBox[{
        SqrtBox["x4"], " ", "x6"}]]}], ",",
     RowBox[{"c0", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SqrtBox["x2"], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3", " ", "x4"}],
       RowBox[{"x1", " ", "x5", " ", "x6"}]]}], ",",
     RowBox[{"m", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"A", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}], ",",
     RowBox[{"mol", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}]}], "}"}]}], "}"}]], "Output",\

 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.55518746845043*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x4"], " ",
         SqrtBox["x5"]}], "x6"]}]}], ",",
     RowBox[{"c0", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SqrtBox["x2"], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{"x4", " ", "x6"}], "x5"]}], ",",
     RowBox[{"m", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ", "x6"}],
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x4"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kg", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}],
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}]]}], ",",
     RowBox[{"s", "\[Rule]",
      FractionBox["1", "x4"]}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x4"], " ",
        SqrtBox["x5"]}], "x6"]}], ",",
     RowBox[{"c0", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SqrtBox["x2"], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{"x4", " ", "x6"}], "x5"]}], ",",
     RowBox[{"m", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ", "x6"}],
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x4"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kg", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}],
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}]]}], ",",
     RowBox[{"s", "\[Rule]",
      FractionBox["1", "x4"]}]}], "}"}]}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468590424*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}],
        RowBox[{
         SqrtBox["x4"], " ", "x6"}]]}]}], ",",
     RowBox[{"c0", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SqrtBox["x2"], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x4", "2"], " ", "x6"}], "x5"]}], ",",
     RowBox[{"m", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"A", "\[Rule]",
      FractionBox["1", "x4"]}], ",",
     RowBox[{"kg", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}],
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}]]}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}],
       RowBox[{
        SqrtBox["x4"], " ", "x6"}]]}], ",",
     RowBox[{"c0", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SqrtBox["x2"], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x4", "2"], " ", "x6"}], "x5"]}], ",",
     RowBox[{"m", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"A", "\[Rule]",
      FractionBox["1", "x4"]}], ",",
     RowBox[{"kg", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}],
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}]]}]}], "}"}]}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468595479*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
     RowBox[{"c0", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SqrtBox["x2"], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{"x4", " ", "x6"}], "x5"]}], ",",
     RowBox[{"kg", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"9", "/", "2"}]], " ",
         SuperscriptBox["x3",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x5"]}],
        RowBox[{
         SuperscriptBox["x1",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x4"], " ",
         SuperscriptBox["x6", "3"]}]]}]}], ",",
     RowBox[{"mol", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x4"], " ",
         SqrtBox["x5"]}],
        RowBox[{
         SqrtBox["x1"], " ", "x6"}]]}]}], ",",
     RowBox[{"s", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}],
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}]]}]}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
     RowBox[{"c0", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SqrtBox["x2"], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{"x4", " ", "x6"}], "x5"]}], ",",
     RowBox[{"kg", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"9", "/", "2"}]], " ",
        SuperscriptBox["x3",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x5"]}],
       RowBox[{
        SuperscriptBox["x1",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x4"], " ",
        SuperscriptBox["x6", "3"]}]]}], ",",
     RowBox[{"mol", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x4"], " ",
        SqrtBox["x5"]}],
       RowBox[{
        SqrtBox["x1"], " ", "x6"}]]}], ",",
     RowBox[{"s", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}],
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}]]}]}], "}"}]}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468600914*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{"{",
   RowBox[{
    RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
    RowBox[{"c0", "\[Rule]", "x2"}], ",",
    RowBox[{"kvK", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3"}],
      RowBox[{"x1", " ", "x6"}]]}], ",",
    RowBox[{"kelvin", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3", " ", "x5"}],
      RowBox[{"x1", " ", "x4", " ",
       SuperscriptBox["x6", "2"]}]]}], ",",
    RowBox[{"kg", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3"}],
      RowBox[{"x1", " ",
       SuperscriptBox["x6", "2"]}]]}], ",",
    RowBox[{"mol", "\[Rule]", "x4"}]}], "}"}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468750674*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{"{",
   RowBox[{
    RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
    RowBox[{"c0", "\[Rule]",
     RowBox[{"x2", " ", "x4"}]}], ",",
    RowBox[{"kvK", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3", " ",
       SuperscriptBox["x4", "2"]}],
      RowBox[{"x1", " ", "x6"}]]}], ",",
    RowBox[{"kelvin", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3", " ", "x4", " ", "x5"}],
      RowBox[{"x1", " ",
       SuperscriptBox["x6", "2"]}]]}], ",",
    RowBox[{"kg", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3"}],
      RowBox[{"x1", " ", "x4", " ",
       SuperscriptBox["x6", "2"]}]]}], ",",
    RowBox[{"s", "\[Rule]",
     FractionBox["1", "x4"]}]}], "}"}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468753119*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
     RowBox[{"c0", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SqrtBox["x2"], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kvK", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SuperscriptBox["x4",
          RowBox[{"3", "/", "2"}]]}],
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"A", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ", "x6"}],
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x4"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kg", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}],
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}]]}], ",",
     RowBox[{"s", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}],
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}]]}]}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
     RowBox[{"c0", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SqrtBox["x2"], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SuperscriptBox["x4",
         RowBox[{"3", "/", "2"}]]}],
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"A", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ", "x6"}],
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x4"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kg", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}],
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}]]}], ",",
     RowBox[{"s", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}],
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}]]}]}], "}"}]}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468755687*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{"{",
   RowBox[{
    RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
    RowBox[{"c0", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "4"], " ", "x3"}],
      RowBox[{"x1", " ",
       SuperscriptBox["x6", "2"]}]]}], ",",
    RowBox[{"kvK", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "9"], " ",
       SuperscriptBox["x3", "3"]}],
      RowBox[{
       SuperscriptBox["x1", "3"], " ",
       SuperscriptBox["x6", "5"]}]]}], ",",
    RowBox[{"kelvin", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "9"], " ",
       SuperscriptBox["x3", "3"], " ", "x5"}],
      RowBox[{
       SuperscriptBox["x1", "3"], " ", "x4", " ",
       SuperscriptBox["x6", "6"]}]]}], ",",
    RowBox[{"mol", "\[Rule]",
     FractionBox[
      RowBox[{"x1", " ", "x4", " ",
       SuperscriptBox["x6", "2"]}],
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3"}]]}], ",",
    RowBox[{"s", "\[Rule]",
     FractionBox[
      RowBox[{"x1", " ",
       SuperscriptBox["x6", "2"]}],
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3"}]]}]}], "}"}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468760353*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
     RowBox[{"c0", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}],
        RowBox[{
         SqrtBox["x2"], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"kvK", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"9", "/", "2"}]], " ",
         SuperscriptBox["x3",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x4"]}],
        RowBox[{
         SuperscriptBox["x1",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x5"], " ",
         SuperscriptBox["x6", "2"]}]]}]}], ",",
     RowBox[{"A", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x1",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x4"], " ",
         SuperscriptBox["x6", "3"]}],
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"9", "/", "2"}]], " ",
         SuperscriptBox["x3",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x5"]}]]}]}], ",",
     RowBox[{"mol", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}], ",",
     RowBox[{"s", "\[Rule]",
      RowBox[{"-",
       FractionBox[
        RowBox[{
         SuperscriptBox["x2",
          RowBox[{"3", "/", "2"}]], " ",
         SqrtBox["x3"], " ",
         SqrtBox["x5"]}],
        RowBox[{
         SqrtBox["x1"], " ",
         SqrtBox["x4"], " ", "x6"}]]}]}]}], "}"}], ",",
   RowBox[{"{",
    RowBox[{
     RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
     RowBox[{"c0", "\[Rule]",
      FractionBox[
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}],
       RowBox[{
        SqrtBox["x2"], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"kvK", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"9", "/", "2"}]], " ",
        SuperscriptBox["x3",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x4"]}],
       RowBox[{
        SuperscriptBox["x1",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x5"], " ",
        SuperscriptBox["x6", "2"]}]]}], ",",
     RowBox[{"A", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x1",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x4"], " ",
        SuperscriptBox["x6", "3"]}],
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"9", "/", "2"}]], " ",
        SuperscriptBox["x3",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x5"]}]]}], ",",
     RowBox[{"mol", "\[Rule]",
      FractionBox[
       RowBox[{"x1", " ", "x4", " ",
        SuperscriptBox["x6", "2"]}],
       RowBox[{
        SuperscriptBox["x2", "3"], " ", "x3"}]]}], ",",
     RowBox[{"s", "\[Rule]",
      FractionBox[
       RowBox[{
        SuperscriptBox["x2",
         RowBox[{"3", "/", "2"}]], " ",
        SqrtBox["x3"], " ",
        SqrtBox["x5"]}],
       RowBox[{
        SqrtBox["x1"], " ",
        SqrtBox["x4"], " ", "x6"}]]}]}], "}"}]}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468895315*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{"{",
   RowBox[{
    RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
    RowBox[{"c0", "\[Rule]", "x2"}], ",",
    RowBox[{"kvK", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "6"], " ",
       SuperscriptBox["x3", "2"]}],
      RowBox[{
       SuperscriptBox["x1", "2"], " ",
       SuperscriptBox["x6", "3"]}]]}], ",",
    RowBox[{"kelvin", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3", " ", "x5"}],
      RowBox[{"x1", " ", "x4", " ",
       SuperscriptBox["x6", "2"]}]]}], ",",
    RowBox[{"A", "\[Rule]",
     FractionBox[
      RowBox[{"x1", " ",
       SuperscriptBox["x6", "2"]}],
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3"}]]}], ",",
    RowBox[{"mol", "\[Rule]",
     FractionBox[
      RowBox[{"x1", " ", "x4", " ",
       SuperscriptBox["x6", "2"]}],
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3"}]]}]}], "}"}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468899894*^9}}],

Cell[BoxData[
 RowBox[{"{",
  RowBox[{"{",
   RowBox[{
    RowBox[{"Rinf", "\[Rule]", "x1"}], ",",
    RowBox[{"c0", "\[Rule]", "x2"}], ",",
    RowBox[{"kvK", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3", " ", "x4"}],
      RowBox[{"x1", " ", "x6"}]]}], ",",
    RowBox[{"kelvin", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3", " ", "x5"}],
      RowBox[{"x1", " ", "x4", " ",
       SuperscriptBox["x6", "2"]}]]}], ",",
    RowBox[{"A", "\[Rule]",
     FractionBox["1", "x4"]}], ",",
    RowBox[{"kg", "\[Rule]",
     FractionBox[
      RowBox[{
       SuperscriptBox["x2", "3"], " ", "x3"}],
      RowBox[{"x1", " ", "x4", " ",
       SuperscriptBox["x6", "2"]}]]}]}], "}"}], "}"}]], "Output",
 CellChangeTimes->{
  3.555186744671112*^9, {3.555186825718958*^9, 3.555186841562527*^9},
   3.555186879293422*^9, 3.555187254679415*^9, 3.555187425311976*^9, {
   3.555187462936023*^9, 3.555187468902075*^9}}]
}, Open  ]]
},
WindowSize->{1212, 844},
WindowMargins->{{Automatic, 153}, {Automatic, 0}},
FrontEndVersion->"8.0 for Linux x86 (32-bit) (October 10, 2011)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[400, 13, 1768, 34, 107, "Input"],
Cell[CellGroupData[{
Cell[2193, 51, 2661, 74, 69, "Input"],
Cell[4857, 127, 683, 12, 33, "Output"],
Cell[5543, 141, 641, 11, 33, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6221, 157, 4061, 110, 183, "Input"],
Cell[10285, 269, 523, 11, 42, "Message"],
Cell[10811, 282, 615, 17, 33, "Output"],
Cell[11429, 301, 523, 11, 42, "Message"],
Cell[11955, 314, 601, 16, 33, "Output"],
Cell[12559, 332, 523, 11, 42, "Message"],
Cell[13085, 345, 574, 15, 30, "Output"],
Cell[13662, 362, 523, 11, 42, "Message"],
Cell[14188, 375, 399, 10, 33, "Output"],
Cell[14590, 387, 523, 11, 42, "Message"],
Cell[15116, 400, 383, 9, 33, "Output"],
Cell[15502, 411, 522, 11, 42, "Message"],
Cell[16027, 424, 613, 17, 33, "Output"],
Cell[16643, 443, 523, 11, 42, "Message"],
Cell[17169, 456, 599, 16, 30, "Output"],
Cell[17771, 474, 398, 10, 33, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[18206, 489, 4530, 134, 126, "Input"],
Cell[22739, 625, 1154, 25, 50, "Output"],
Cell[23896, 652, 1024, 21, 50, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24957, 678, 5176, 152, 259, "Input"],
Cell[30136, 832, 2710, 91, 108, "Output"],
Cell[32849, 925, 3397, 112, 108, "Output"],
Cell[36249, 1039, 3181, 106, 108, "Output"],
Cell[39433, 1147, 2741, 91, 108, "Output"],
Cell[42177, 1240, 2838, 95, 108, "Output"],
Cell[45018, 1337, 3160, 106, 108, "Output"],
Cell[48181, 1445, 964, 27, 50, "Output"],
Cell[49148, 1474, 1064, 30, 50, "Output"],
Cell[50215, 1506, 3188, 108, 108, "Output"],
Cell[53406, 1616, 1397, 42, 50, "Output"],
Cell[54806, 1660, 3589, 118, 108, "Output"],
Cell[58398, 1780, 1200, 35, 50, "Output"],
Cell[59601, 1817, 1008, 28, 50, "Output"]
}, Open  ]]
}
]
*)
