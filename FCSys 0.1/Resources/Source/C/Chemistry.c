// Functions to read information from chemical formulas
//
// See the corresponding functions in the FCSys.BaseClasses.Utilities package
// (Modelica) for full descriptions and examples.
//
// Copyright (C) 2012, Kevin Davies.
//
// The content of this file is free software; it can be redistributed and/or
// modified under the terms of the Modelica License 2, see the license
// conditions and the accompanying disclaimer in file
// FCSys/resources/documentation/ModelicaLicense2.html or in
// FCSys.UsersGuide.ModelicaLicense2.
// Initial version: 2012-10-07

// Prevent multiple inclusion.
#ifndef _FCSYS_CHEMISTRY_INCLUDED_
#define _FCSYS_CHEMISTRY_INCLUDED_

#include <ctype.h>
#include <string.h>
#include "ModelicaUtilities.h"

static int charge(const char* formula)
{
    // Return the charge of a species based on its chemical formula.

    int z = 0; // Charge number
    int i = 0; // Index

    while (formula[i] != '\0') {

        // Ignore leading whitespace.
        while (formula[i] != '\0' && isspace(formula[i]))
            i++;

        // Check that the symbol begins with a letter.
        if (!isalpha(formula[i])) {
            return 0; // Error
            // Currently, an error cannot be distinguished from zero charge.
        }

        // Advance past the symbol.  It may continue with lowercase letters.
        while (formula[++i] != '\0' && islower(formula[i]));

        // Advance past the coefficient.
        while (formula[i] != '\0' && isdigit(formula[i]))
            i++;

        // Read the charge.
        if (formula[i] != '\0' && (formula[i] == '+' || formula[i] == '-')) {
            if (isdigit(formula[i+1])) {
                z += atoi(formula + i); // Continues until nondigit
                while (formula[++i] != '\0' && isdigit(formula[i]));
            } else if (formula[i] == '+') {
                z++;
                i++;
            } else if (formula[i++] == '-')
                z--;
        }
    }

    return z;
}

static int countElements(const char* formula)
{
    // Return the number of elements in a species based on its chemical formula.

    int z = 0; // Charge number
    int i = 0; // Index
    int n = 0; // Number of unique elements

    while (formula[i] != '\0') {

        // Ignore leading whitespace.
        while (formula[i] != '\0' && isspace(formula[i]))
            i++;

        // Interpret the symbol.
        if (formula[i] == 'e') {
            n--; // Electrons are counted via their charge (below).
            i++;
        } else if (isalpha(formula[i])) {
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
            if (isdigit(formula[i+1])) {
                z += atoi(formula + i); // Continues until nondigit
                while (formula[++i] != '\0' && isdigit(formula[i]));
            } else if (formula[i] == '+') {
                z++;
                i++;
            } else if (formula[i++] == '-')
                z--;
        }
        n++;
    }

    // Count electrons as a present element (or rather particle) if the charge
    // is nonzero.
    if (z != 0)
        n++;

    return n;
}

static void readElement(const char* formula, char** symbol, int* n, int* z,
                        char** remainder)
{
    // Read the symbol, coefficient, and charge of an element from a chemical
    // formula.  Also return the remainder of the formula

    int i, symbol_length;

    // Ignore leading whitespace.
    int symbol_start = ModelicaStrings_skipWhiteSpace(formula, 1) - 1;

    // Check that the symbol begins with a letter.
    if (!isalpha(formula[symbol_start])) {
        *symbol = ""; // Indicate an error.
        return;
    }

    // The symbol may continue with lowercase letters.
    i = symbol_start + 1; // Working index
    while (formula[i] != '\0' && islower(formula[i]))
        i++;

    // Copy the symbol (return by reference).
    symbol_length = i - symbol_start;
    *symbol = ModelicaAllocateString(symbol_length);
    strncpy(*symbol, &formula[symbol_start], symbol_length);

    // Read the coefficient.
    if (formula[i] != '\0' && isdigit(formula[i])) {
        *n = atoi(formula + i); // Continues until nondigit
        while (formula[++i] != '\0' && isdigit(formula[i]));
    } else
        *n = 1;

    // Read the charge.
    if (formula[i] != '\0' && (formula[i] == '+' || formula[i] == '-')) {
        if (isdigit(formula[i+1])) {
            *z = atoi(formula + i); // Continues until nondigit
            while (formula[++i] != '\0' && isdigit(formula[i]));
        } else if (formula[i] == '+') {
            *z = 1;
            i++;
        } else if (formula[i++] == '-')
            *z = -1;
    } else
        *z = 0;

    // Return the remainer too.
    *remainder = &formula[i];
    return;
}

#endif /* _FCSYS_CHEMISTRY_INCLUDED_ */
