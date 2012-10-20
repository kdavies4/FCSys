/* Functions to read information from chemical formulas

   See the corresponding functions in the FCSys.BaseClasses.Utilities package
   (Modelica) for full descriptions and examples.

   Copyright (C) 2012, Kevin Davies.

   The content of this file is free software; it can be redistributed
   and/or modified under the terms of the Modelica License 2, see the
   license conditions and the accompanying disclaimer in file
   FCSys/resources/documentation/ModelicaLicense2.html or in
   FCSys.UsersGuide.ModelicaLicense2.                                         */
// Initial version: 2012-10-07
#include <ctype.h>
#include <string.h>
#include "ModelicaUtilities.h"

static int charge(const char* formula)
{
    /* Return the charge of a species based on its chemical formula.          */

    int charge = 0;
    int i = 0; // Index

    while (formula[i] != '\0') {

        // Ignore leading whitespace.
        while (formula[i] != '\0' && isspace(formula[i]))
            i++;

        // Check that the symbol begins with a letter.
        if (!isalpha(formula[i])) {
            return 0; // Error
            // Currently, an error is not distinguished from zero charge.
        }

        // Advance past the symbol.  It may continue with lowercase letters.
        while (formula[++i] != '\0' && islower(formula[i]));

        // Advance past the coefficient.
        while (formula[i] != '\0' && isdigit(formula[i]))
            i++;

        // Read the charge.
        if (formula[i] != '\0' && (formula[i] == '+' || formula[i] == '-')) {
            if isdigit(formula[i+1]){
                charge += atoi(formula + i); // automatically continues until nondigit
                while (formula[++i] != '\0' && isdigit(formula[i]));
            } else if (formula[i] == '+'){
                charge++;
                i++;
            } else if (formula[i++] == '-')
                charge--;
        }
    }

    return charge;
}

static int countElements(const char* formula)
{
    /* Return the number of unique elements in a species based on its chemical
       formula.                                                               */

    int charge = 0;
    int i = 0; // Index
    int n = 0; // Number of unique elements

    while (formula[i] != '\0') {

        // Ignore leading whitespace.
        while (formula[i] != '\0' && isspace(formula[i]))
            i++;

        // Interpret the symbol.
        if (formula[i] == 'e'){
            n--; // Electrons are counted via their charge (below).
            i++;
        } else if isalpha(formula[i]) {
            // Advance past the symbol.  It may continue with lowercase letters.
            while (formula[++i] != '\0' && islower(formula[i]));
        } else {
            return 0; // Error; the symbol does not start with a letter.
            // Currently, an error is not distinguished from zero elements.
        }

        // Advance past the coefficient.
        while (formula[i] != '\0' && isdigit(formula[i]))
            i++;

        // Read the charge.
        if (formula[i] != '\0' && (formula[i] == '+' || formula[i] == '-')) {
            if isdigit(formula[i+1]){
                charge += atoi(formula + i); // automatically continues until nondigit
                while (formula[++i] != '\0' && isdigit(formula[i]));
            } else if (formula[i] == '+'){
                charge++;
                i++;
            } else if (formula[i++] == '-')
                charge--;
        }
        n++;
    }

    // Count electrons as a present element (or rather particle) if the charge
    // is nonzero.
    if (charge != 0)
        n++;

    return n;
}

static void readElement(const char* formula, int startIndex, char** symbol,
                        int* coeff, int* charge, int* nextIndex)
{
    /* Read the symbol, coefficient, and charge of an element at an index in a
       chemical formula.

       Note: startIndex and nextIndex are 1-based indices (Modelica format), but
       all other indices are 0-based (C format).                              */

    // Ignore leading whitespace.
    int symbol_start = ModelicaStrings_skipWhiteSpace(formula, startIndex) - 1;

    // Check that the symbol begins with a letter.
    if (!isalpha(formula[symbol_start])) {
        *nextIndex = 0; // Error
        return;
    }

    // The symbol may continue with lowercase letters.
    int i = symbol_start + 1; // Working index
    while (formula[i] != '\0' && islower(formula[i]))
        i++;

    // Copy the symbol (return by reference).
    int symbol_length = i - symbol_start;
    *symbol = ModelicaAllocateString(symbol_length);
    strncpy(*symbol, &formula[symbol_start], symbol_length);

    // Read the coefficient.
    if (formula[i] != '\0' && isdigit(formula[i])) {
        *coeff = atoi(formula + i); // automatically continues until nondigit
        while (formula[++i] != '\0' && isdigit(formula[i]));
    } else
        *coeff = 1;

    // Read the charge.
    if (formula[i] != '\0' && (formula[i] == '+' || formula[i] == '-')) {
        if isdigit(formula[i+1]){
            *charge = atoi(formula + i); // automatically continues until nondigit
            while (formula[++i] != '\0' && isdigit(formula[i]));
        } else if (formula[i] == '+'){
            *charge = 1;
            i++;
        } else if (formula[i++] == '-')
            *charge = -1;
    } else
        *charge = 0;

    // Return the index by reference.
    *nextIndex = i + 1; // plus one for Modelica format
    return;
}
