(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 8.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc.                                               *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[       835,         17]
NotebookDataLength[      4173,        126]
NotebookOptionsPosition[      4505,        119]
NotebookOutlinePosition[      4839,        134]
CellTagsIndexPosition[      4796,        131]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ",
   RowBox[{
   "Calculate", " ", "the", " ", "transcendental", " ", "numbers", " ", "in",
    " ", "the", " ", "Wien", " ", "displacement", " ", "law", " ",
    "constants", " ",
    RowBox[{"(",
     RowBox[{"see", " ",
      RowBox[{"Units", ".", "mo"}]}], ")"}]}], " ", "*)"}],
  "\[IndentingNewLine]",
  RowBox[{
   RowBox[{
    RowBox[{"c3", "[", "u_", "]"}], ":=",
    RowBox[{"NSolve", "[",
     RowBox[{
      RowBox[{
       RowBox[{
        RowBox[{
         RowBox[{"Exp", "[", "x", "]"}], "*",
         RowBox[{"(",
          RowBox[{"u", "-", "x"}], ")"}]}], "\[Equal]", "u"}], "&&",
       RowBox[{"x", ">", "0"}]}], ",", "x"}], "]"}]}], "\[IndentingNewLine]",

   RowBox[{"(*", " ",
    RowBox[{
    "Wien", " ", "frequency", " ", "displacement", " ", "law", " ",
     "constant"}], " ", "*)"}], "\[IndentingNewLine]",
   RowBox[{"NumberForm", "[",
    RowBox[{
     RowBox[{"c3", "[", "3", "]"}], ",",
     RowBox[{"{",
      RowBox[{"16", ",", "15"}], "}"}]}], "]"}], "\[IndentingNewLine]",
   RowBox[{"(*", " ",
    RowBox[{
    "For", " ", "the", " ", "calculation", " ", "of", " ", "the", " ", "Wien",
      " ", "wavelength", " ", "displacement", " ", "law", " ", "constant",
     " ",
     RowBox[{"(",
      RowBox[{"see", " ",
       RowBox[{"Units", ".", "mo"}]}], ")"}]}], " ", "*)"}],
   "\[IndentingNewLine]",
   RowBox[{"NumberForm", "[",
    RowBox[{
     RowBox[{"c3", "[", "5", "]"}], ",",
     RowBox[{"{",
      RowBox[{"16", ",", "15"}], "}"}]}], "]"}]}]}]], "Input",
 CellChangeTimes->{{3.547804858810236*^9, 3.547804903847475*^9}, {
  3.547805195959431*^9, 3.547805206123515*^9}, {3.547805430464469*^9,
  3.547805541036757*^9}, {3.547805586255049*^9, 3.547805633025752*^9}, {
  3.547805676115864*^9, 3.547806092346467*^9}, {3.54780614099026*^9,
  3.547806214564818*^9}, {3.547816808138461*^9, 3.54781692408461*^9}, {
  3.547816954692015*^9, 3.547816955150472*^9}, {3.547830644565884*^9,
  3.547830657949221*^9}, {3.547830704961005*^9, 3.547830722229373*^9}, {
  3.547830812742705*^9, 3.547830824358557*^9}}],

Cell[BoxData[
 TagBox[
  RowBox[{"{",
   RowBox[{"{",
    RowBox[{"x", "\[Rule]",
     InterpretationBox["\<\"2.821439372122079\"\>",
      2.8214393721220787`,
      AutoDelete->True]}], "}"}], "}"}],
  NumberForm[#, {16, 15}]& ]], "Output",
 CellChangeTimes->{
  3.547805975614977*^9, {3.547806014757172*^9, 3.5478060940034*^9}, {
   3.547806159400521*^9, 3.547806219863239*^9}, {3.547816823614075*^9,
   3.547816924692675*^9}, {3.547830650663012*^9, 3.547830660088499*^9},
   3.547830723890523*^9, {3.547830818328186*^9, 3.547830825825211*^9}}],

Cell[BoxData[
 TagBox[
  RowBox[{"{",
   RowBox[{"{",
    RowBox[{"x", "\[Rule]",
     InterpretationBox["\<\"4.965114231744276\"\>",
      4.965114231744276,
      AutoDelete->True]}], "}"}], "}"}],
  NumberForm[#, {16, 15}]& ]], "Output",
 CellChangeTimes->{
  3.547805975614977*^9, {3.547806014757172*^9, 3.5478060940034*^9}, {
   3.547806159400521*^9, 3.547806219863239*^9}, {3.547816823614075*^9,
   3.547816924692675*^9}, {3.547830650663012*^9, 3.547830660088499*^9},
   3.547830723890523*^9, {3.547830818328186*^9, 3.5478308258768*^9}}]
}, Open  ]]
},
WindowSize->{740, 857},
WindowMargins->{{0, Automatic}, {Automatic, 0}},
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
Cell[CellGroupData[{
Cell[1257, 32, 2126, 54, 164, "Input"],
Cell[3386, 88, 552, 13, 43, "Output"],
Cell[3941, 103, 548, 13, 43, "Output"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature jvT8osyNM2Q#gDwgrmU1bcDG *)
