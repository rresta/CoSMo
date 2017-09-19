#include <iostream>
#include <fstream>
#include <ostream>
#include <string>

#include <cmath>
#include <vector>
#include <sstream>
#include <cstdio>
clock_t start = clock();

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

using namespace seqan;

int main () {
    std::string line;

    //tools that you want to use

    std::string tool_name_1("RNAfold");
    std::string tool_name_2("RME");
    std::string tool_name_3("RNAstructure");
    std::string tool_name_4("SeqFold");

    std::ifstream myfile("../data/config_file.txt");
    if (!(myfile.is_open()))       // reading a configuration file
    {
        std::cout << "Unable to open configuration file";
        exit (EXIT_FAILURE);
    }

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
        // I modified the file name in order to have the dbn extention at the end.
        // The original extention of the output tool was .fold .


        RnaStructFileIn inputfile;
        if (!open(inputfile, "../data/in.dbn", OPEN_RDONLY))            //This inputfile is the output of the first tool
            // if (!open(inputfile, "1out_Human_tRNA_38.chr6.fold", OPEN_RDONLY))
        {
            std::cerr << "could not open inputfile" << std::endl;
            return 1;
        }
        RnaStructFileIn inputfile2;                             //This inputfile is the output of the second tool
        if (!open(inputfile2, "../data/in2.dbn", OPEN_RDONLY))
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

        //one record of RnaStructContents is a sequence
        RnaStructContents contents{};
        RnaStructContents contents2{};
        readRecords(contents, inputfile, 100);
        readRecords(contents2, inputfile2, 100);

        //////////////////////////////////////////////////////////////////////////////////////////////
        //Compare sequences
        //////////////////////////////////////////////////////////////////////////////////////////////
        //Compare the sequnces of the files obtained by using the different tools
        //If the sequence is the same, the sequence is written in an output file.

        //NEXT: crea un nuovo grafo(stringa--) e appendici il primo e il secondo.

        if (contents.records[0].sequence== contents2.records[0].sequence) {
            writeRecords(outputfile, contents);          //AFTER -APPEND FUNCTION- THE OUTPUT FILE IS EMPTY.NO WRITTEN RECORDS.NO ERRORS
            std::cout << "The sequence of the files obtained by the tools is the same." << std::endl;
            std::cout << "Sequence: " <<contents.records[0].sequence << std::endl;
            std::cout << "Length: " <<contents.records[0].seqLen << std::endl;
            //std::cout << "Graph in 0 position of the FixedGraph String: " << contents.records[0].fixedGraphs[0];
        }
        else {
            std::cout << "The sequence if the files is not the same." << std::endl;
            exit (EXIT_FAILURE);
        }

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
                std::cout << vett[m] << " ";
        }
        //std::cout << "Dimension of the vector that contains the graph of the first file: " << length(vett);

        RnaAdjacencyIterator adj_it(contents.records[0].fixedGraphs[0].inter, 0);
        std::cout << "pos 0: " << " " << value(adj_it) << std::endl;  //***if you check on the output.ebpseq this is not the 0 but the 1 indexing positions

        //////////////////////////////////////////////////////////////////////////////////////////////
        //Creation of the matrix with the base-pairs and the respoctive frequences
        //////////////////////////////////////////////////////////////////////////////////////////////

        //These will be the arrays obtained in the previous part
//        int RNAfold_vett[contents.records[0].seqLen];            /////////////////////////////////////////////////////////////////////
//        int RME_vett[contents.records[0].seqLen];                /////////////////////////////////////////////////////////////////////
//        int RNAstructure_vett[contents.records[0].seqLen];       /////////////////////////////////////////////////////////////////////
//        int SeqFold_vett[contents.records[0].seqLen];            /////////////////////////////////////////////////////////////////////

        //matrix[n][m]: n->sequnce length m->(number of tools=4) x 2-> in order to have for each possible value the frequence
//        int matrix[contents.records[0].seqLen][8];

        int RME_flag=0;
        int RNAstructure_flag=0;
        int SeqFold_flag=0;

        int RNAfold_count=1;
        int RME_count=1;
        int RNAstructure_count=1;
        int SeqFold_count=1;

        /////////////////////////////////////////////////////////////////////
        //Examples of arrays obtained by accessing the graphs
        //The length of the sequence in this case is 3.
        /////////////////////////////////////////////////////////////////////
        int RNAfold_vett[3]={5, 6, 7};
        int RME_vett[3]={9,6,7};
        int RNAstructure_vett[3]={9, 11, 7};
        int SeqFold_vett[3]={5, 11, 7};
        int matrix[3][8];
        /////////////////////////////////////////////////////////////////////

 //       for (int n=0; n<contents.records[0].seqLen; n++) {                          /////////////////////////////////////////////////////////////////////
        for (int n=0; n<3; n++) {
            int m = 0;
            matrix[n][m] = RNAfold_vett[n];
            if (RNAfold_vett[n] == RME_vett[n]) {
                RNAfold_count++;
                RME_flag = 1;
            } else {
                m = m + 2;
                matrix[n][m] = RME_vett[n];
                if (RME_vett[n]==RNAstructure_vett[n]) {
                    RME_count++;
                    RNAstructure_flag=1;
                }
                if (RME_vett[n]==SeqFold_vett[n]) {
                    RME_count++;
                    SeqFold_flag=1;
                }
                matrix[n][m+1] = RME_count;
            }
            if (RNAfold_vett[n] == RNAstructure_vett[n]) {
                RNAfold_count++;
                RNAstructure_flag = 1;
            } else if (RNAstructure_flag != 1) {
                m = m + 2;
                matrix[n][m] = RNAstructure_vett[n];
                if (RNAstructure_vett[n]==SeqFold_vett[n]) {
                    RNAstructure_count++;
                    SeqFold_flag=1;
                }
                matrix[n][m+1] = RNAstructure_count;
            }
            if (RNAfold_vett[n] == SeqFold_vett[n]) {
                RNAfold_count++;
                SeqFold_flag = 1;
            } else if (SeqFold_flag != 1 ){
                m = m + 2;
                matrix[n][m] = SeqFold_vett[n];
                matrix[n][m+1] = SeqFold_count;
            }
            matrix[n][1] = RNAfold_count; ///////////////////////////////////////

            std::cout << "> " << n+1 << " nucleotide : " << std::endl;
            std::cout << "RNAfold COUNT: " << RNAfold_count << ",  ";
            std::cout << "RME COUNT: " << RME_count << ",  ";
            std::cout << "RNAstructure COUNT: " << RNAstructure_count << ",  ";
            std::cout << "SeqFold COUNT: " << SeqFold_count << std::endl;

            RNAfold_count = 1;
            RME_count = 1;
            RNAstructure_count = 1;
            SeqFold_count = 1;

            RME_flag = 0;
            RNAstructure_flag = 0;
            SeqFold_flag = 0;
        }

        std::cout << "MATRIX:" << std::endl;
        for (int j=0; j<3; j++){
            for (int k=0; k<8; k++){
                std::cout << matrix[j][k] << " ";
            }
            std::cout << std::endl;
        }
        /////////////////////////////////////////////////////////////////////
        //SHAPE
        /////////////////////////////////////////////////////////////////////
        //SHAPE data contribution to have the right weights

        /////////////////////////////////////////////////////////////////////
        //Graph
        /////////////////////////////////////////////////////////////////////

    //Example:
//    RnaStructureGraph bgraph;
//    bgraph.specs = "The bpp structure";
//    for (typename Size<Rna5String>::Type idx = 0u; idx < contents.records[0].seqLen; ++idx)
//        addVertex(bgraph.inter);
//    addEdge(bgraph.inter, 0u, 2u, 0.8);
//    addEdge(bgraph.inter, 1u, 2u, 0.2);
//    append(contents.records[0].bppMatrGraphs, bgraph);

    myfile.close();

    clock_t end = clock();
    double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;

    std::cout << "\nTIME (sec): " << time;

    return 0;
}