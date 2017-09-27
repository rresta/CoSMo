#include <iostream>
#include <fstream>
#include <ostream>
#include <string>

#include <cmath>
#include <vector>
#include <sstream>
#include <cstdio>
#include <limits>
clock_t start = clock();

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/arg_parse.h>
#include <seqan/store.h>
#include <seqan/stream.h>
#include <seqan/align.h>
#include <seqan/align_rna.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/rna_io.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

// defines all the constants used in the app
//#include "data_types.h"
#include "option.h"

using namespace seqan;

int main (int argc, char const ** argv)
    {
        // Argument Parser
        if (argc == 1)
        {
            std::cout << "type " << argv[0] << " --help to get the parameters table (-i option is mandatory)" << std::endl;
            return 1;
        }

        ArgumentParser parser;
        Options options;
        setupArgumentParser(parser, options);
        ArgumentParser::ParseResult res;
        res = parse(options, parser, argc, argv); // Fill the options structure
        if (res != ArgumentParser::ParseResult::PARSE_OK)
            return res == ArgumentParser::ParseResult::PARSE_ERROR ? 1 : 0;

        // Read input files

        RnaStructFileIn inputfile;                             //This inputfile is the output of the second tool
        if (!open(inputfile, toCString(options.inFile), OPEN_RDONLY))
            // if (!open(inputfile, "1out_Human_tRNA_38.chr6.fold", OPEN_RDONLY))
        {
            std::cerr << "could not open inputfile" << std::endl;
            return 1;
        }

        RnaStructFileOut outputfile;
        if (!open(outputfile, "../data/output.ebpseq", OPEN_RDWR)) {
            std::cerr << "could not open outputfile" << std::endl;
            return 1;
        }

        RnaStructContents contents;
        readRecords(contents, inputfile, std::numeric_limits<unsigned>::max());
        _VVV(options, contents.header.description);
        for(unsigned i = 0; i < length(contents); ++i)
        {
            _VVV(options, "iteration i = " << i);
            _VVV(options, contents.records[i].sequence);
            _VVV(options,contents.records[i].name);
            for(unsigned j = 0; j < length(contents.records[i].fixedGraphs); ++j)
            {
                _VVV(options, "iteration j = " << j);
                _VVV(options, contents.records[i].fixedGraphs[j].inter);
                _VVV(options, contents.records[i].fixedGraphs[j].specs);
                _VVV(options, contents.records[i].fixedGraphs[j].energy);
            }
        }

        return 167;
//        _readRnaInputFile(filecontents1, options.inFile, options);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //Run the line of the configuration file
    //////////////////////////////////////////////////////////////////////////////////////////////

//        while ( getline (myfile,line) )
//        {
//            cout << '\n' << line;

        //check if the names of the tools are in the lines and run RNAfold

//            size_t found1 = line.find(tool_name_1);      //Unsigned integral type
//            if (found1!=string::npos)
//            {
//                cout << "\nThis row is RNAfold\n";
//                int i;
//                printf ("Checking if processor is available...");
//                if (system(NULL)) puts ("Ok");
//                else exit (EXIT_FAILURE);
//                printf ("Executing command: RNAfold --outfile=1out < trna.fasta\n");
//                i=system ("RNAfold --outfile=1out < trna.fasta");  //the file has to be in the /home/rossellar/CLionProjects/prova_funzioni/cmake-build-debug
//            }
//        }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //Open files
    //////////////////////////////////////////////////////////////////////////////////////////////
    //the input file has been obtained by using RNAfold on Ubuntu system.

//    char oldname[] ="../data/1out_Human_tRNA_2.chr7.fold";
//    char newname[] ="../data/in.dbn";
//    int result =rename( oldname , newname );
//    if (result!= 0)
//    {
//        std::cerr << "ERROR: the name of the output file of RNAfold is not valid." << std::endl;
//        return 1;
//    }


    //////////////////////////////////////////////////////////////////////////////////////////////
    //Compare sequences
    //////////////////////////////////////////////////////////////////////////////////////////////
    //Compare the sequnces of the files obtained by using the different tools
    //If the sequence is the same, the sequence is written in an output file.
    RnaStructContents contents2;

    if (contents.records[0].sequence== contents2.records[0].sequence) {
        writeRecords(outputfile, contents);          //AFTER -APPEND FUNCTION- THE OUTPUT FILE IS EMPTY.NO WRITTEN RECORDS.NO ERRORS
        std::cout << "The sequence of the files obtained by the tools is the same." << std::endl;
        std::cout << "Sequence: " <<contents.records[0].sequence << std::endl;
        std::cout << "Length: " <<contents.records[0].seqLen << std::endl;
        std::cout << "-> Sequence position access: " << contents.records[0].sequence[0] << std::endl;
        //std::cout << "Graph in 0 position of the FixedGraph String: " << contents.records[0].fixedGraphs[0];
    }
    else {
        std::cout << "The sequence if the files is not the same." << std::endl;
        exit (EXIT_FAILURE);
    }
    return 0;
    //////////////////////////////////////////////////////////////////////////////////////////////
    //Definition of null objects
    //////////////////////////////////////////////////////////////////////////////////////////////
    //MISSING: count the number of neighbors, if they are 0 the corresponding base-pair will have 0.

    //////////////////////////////////////////////////////////////////////////////////////////////
    //Access to the graphs obtained by using the tools
    //////////////////////////////////////////////////////////////////////////////////////////////
    //Access to the graphs by using 'RnaAdjacencyIterator' and 'value' from SeqAn library and
    //creation of the arrays that will cointain the graphs.
    //PROBLEM: shifted positions***

//      for (int m=0; m<contents.records[0].seqLen; m++) {      //  -> does not work with RnaAdjacencyIterator, it needs the check on neighbors before.

    int vett[7];
    for (int m=0; m<7; m++) {                             //It works for the first 8 positions
            RnaAdjacencyIterator adj_it(contents.records[0].fixedGraphs[0].inter, m);
            vett[m]=value(adj_it);
            //std::cout << vett[m] << " ";
    }

    std::cout << std::endl << contents.records[0].fixedGraphs[0].inter;

    /////////////////////////////////////////////////////////////////////
    //OutputGraph creation by copying the first one
    /////////////////////////////////////////////////////////////////////
/*
    typedef unsigned int TCargo;
    typedef Graph<Undirected<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    TGraph OutputGraph;

    for (int i=0; i<contents.records[0].seqLen; i++)    {
        TVertexDescriptor contents.records[0].sequence[i] = addVertex(OutputGraph);
    }
    for (int i=0; i<contents.records[0].seqLen; i++)     {
        addEdge(OutputGraph, , , );
    }


    */

    /////////////////////////////////////////////////////////////////////
    //Updating of the output graph by accessing to the others(if the edge exists->updating of the weight, if it does not exists->creation of the edge)
    /////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////
    //SHAPE
    /////////////////////////////////////////////////////////////////////
    //SHAPE data contribution to have the right weights

    /////////////////////////////////////////////////////////////////////
    //Filtering of the bpp matrix by using the output graph
    /////////////////////////////////////////////////////////////////////

//Example:
//    RnaStructureGraph bgraph;
//    bgraph.specs = "The bpp structure";
//    for (typename Size<Rna5String>::Type idx = 0u; idx < contents.records[0].seqLen; ++idx)
//        addVertex(bgraph.inter);
//    addEdge(bgraph.inter, 0u, 2u, 0.8);
//    addEdge(bgraph.inter, 1u, 2u, 0.2);
//    append(contents.records[0].bppMatrGraphs, bgraph);


    clock_t end = clock();
    double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;

    std::cout << "\nTIME (sec): " << time;

    return 0;
}