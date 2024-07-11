/* ************************************************************
 * CS476 Project #1 solution                                  *
 * 							                                  *
 * Your Name: Trevor Wright                                   *
 * Your last-three:  560                                      *
 * Your course section #: 001                                 *
 *                                                            *
 * Summer 2024                                                *
 *                                                            *
 *                                                            *
 * ************************************************************/

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

// Helper functions
string getSequenceFromFile(string filename);                                                                // gather gene sequences from file
vector<vector<float>> getSubstitutionMatrix(string filename, string& subMatrixName);                        // gets a substitution matrix from a file
int determineLocation(char letter);                                                                         // determine location of letter in substitution matrix

// Sequence Alignment functions
vector<vector<float>> globalAlignment(string gene1, string gene2, vector<vector<float>> subMatrix, float gapPenalty, string& optPath);
vector<vector<float>> localAlignment(string gene1, string gene2, vector<vector<float>> subMatrix, float gapPenalty, string& optPath, string& alignedGene1, string& alignedGene2, int& rowMax, int& columnMax);
vector<vector<float>> semiAlignment(string gene1, string gene2, vector<vector<float>> subMatrix, float gapPenalty, string& optPath, string& alignedGene1, string& alignedGene2, int& rowMax, int& columnMax);

/* ****************************************************************
 * steps:                                                         *
 * 1. Initialize OPT and Path Table                               *
 * 2. Use DP to determine Optimal Path                            *
 * 3. During DP insert D, V, H into Path matrix                   *
 * 4. Selecting an index inside of Path matrix traverse in        *
 *    reverse to figure out forward path                          *
 * 5. reverse forward path to have it ouput correctly             * 
 * ****************************************************************/


// Output to console
void display(string& choice1, string& choice2, string& choice3, string& choice4, float& gapPenalty);        // prints interface for user
void printSubstitutionMatrix(vector<vector<float>> arr, string subMatrixName);                              // prints substitution matrix
void printAlignedSequence(string gene1, string gene2, string path, string type);                            // prints the alinged sequence


int main(){
    fstream fin;
    string choice1, choice2, choice3, choice4, gene1, gene2, subMatrixName, optPath, alignedGene1, alignedGene2;
    float gapPenalty;
    int rowMax = 0, columnMax = 0;
    vector<vector<float>> subMatrix;
    vector<vector<float>> OPT;

    display(choice1, choice2, choice3, choice4, gapPenalty);

    cout << "\nfilename for Amino Acid Sequenece 1: " << choice1 << endl;
    cout << "filename for Amino Acid Sequenece 2: " << choice2 << endl;
    cout << "filename for Substitution Matrix: " << choice3 << endl;

    gene1 = getSequenceFromFile(choice1);
    gene2 = getSequenceFromFile(choice2);

    cout << endl;
    cout << "Gene 1: " << gene1 << endl;
    cout << "Gene 2: " << gene2 << endl;
    cout << "Method: " << choice4 << endl;;
    cout << "Gap Penalty: " << gapPenalty << endl;
    cout << endl;

    subMatrix = getSubstitutionMatrix(choice3, subMatrixName);

    cout << endl;

    printSubstitutionMatrix(subMatrix, subMatrixName);

    cout << endl;

    if(choice4 == "global") OPT = globalAlignment(gene1, gene2, subMatrix, gapPenalty, optPath);
    else if (choice4 == "local") OPT = localAlignment(gene1, gene2, subMatrix, gapPenalty, optPath, alignedGene1, alignedGene2, rowMax, columnMax);
    else OPT = semiAlignment(gene1, gene2, subMatrix, gapPenalty, optPath, alignedGene1, alignedGene2, rowMax, columnMax);

    cout << "\t-\t";
    for(int i = 0; i < gene1.length(); i++){
        cout << gene1[i] << "\t";
    }

    cout << endl;

    for (int i = 0; i <= gene2.length(); ++i) {
        if(i == 0) cout << "-\t";

        else cout << gene2[i - 1] << "\t";
        for (int j = 0; j <= gene1.size(); ++j) {
            cout << OPT[i][j] << "\t";
        }
        cout << endl;
    }

    fin.close();

    reverse(optPath.begin(), optPath.end());
    cout << endl;
    if(choice4 == "global") cout << "Alignment Score: " << OPT[gene2.length()][gene1.length()] << endl;
    else if(choice4 == "local") cout << "Alignment Score: " << OPT[rowMax][columnMax] << endl;
    else cout << "Alignment Score: " << OPT[rowMax][columnMax] << endl;
    cout << "Forward Path: " << optPath << endl;
    cout << "Aligned Sequences: " << endl;
    if(choice4 == "global") printAlignedSequence(gene1, gene2, optPath, choice4);
    else if(choice4 == "local") printAlignedSequence(alignedGene1, alignedGene2, optPath, choice4);
    else printAlignedSequence(alignedGene1, alignedGene2, optPath, choice4);

    return 0;
}

void display(string& choice1, string& choice2, string& choice3, string& choice4, float& gapPenalty){
    char option;

    cout << "Enter the name of the first amino acide sequence file: ";
    cin >> choice1;

    cout << endl;
    cout << "Enter the name of the second amino acid sequence file: ";
    cin >> choice2;

    cout << endl;
    cout << "Enter the name of the substitution matrix: ";
    cin >> choice3;

    cout << endl;
    cout << "Select the alignment type" << endl;
    cout << "\t (A) Global" << endl;
    cout << "\t (B) Local" << endl;
    cout << "\t (C) Semi-Global" << endl;

    cout << "Option: ";
    cin >> option;

    switch(toupper(option)){
        case 'A':
            choice4 = "global";
            break;
        case 'B':
            choice4 = "local";
            break;
        case 'C':
            choice4 = "semi";
            break;
        default:
            cout << "selected non of the options. Default is Global" << endl;
            choice4 = "global";
            break;
    }

    cout << endl;
    cout << "Emter a number to set the gap penalty: ";
    cin >> gapPenalty;

}

string getSequenceFromFile(string filename){
    fstream fin;
    string gene;

    fin.open(filename);

    if(!fin.is_open()){
        cerr << "The file is unable to open " << filename << endl;
    }

    getline(fin, gene);
    getline(fin,gene);

    fin.close();

    return gene;
}

vector<vector<float>> getSubstitutionMatrix(string filename, string& subMatrixName){
    fstream fin;
    vector<vector<float>> matrix;
    string skip, name, MatrixName;
    char letter;

    fin.open(filename);

    if(!fin.is_open()){
        cerr << "This file was unable to open " << filename << endl;
    }

    getline(fin, skip);
    getline(fin, skip);

    stringstream ws(skip);
    while(ws >> letter){
        name += letter;
        if(ws.peek() == ',') ws.ignore();
    }

    while(getline(fin, skip)){
        if(skip.empty()) break;

        vector<float> row;
        stringstream ss(skip);
        float value;

        while(ss >> value){
            row.push_back(value);
            if(ss.peek() == ',') ss.ignore();
        }

        matrix.push_back(row);

    }

    fin.close();

    subMatrixName = name;

    return matrix;
}

void printSubstitutionMatrix(vector<vector<float>> arr, string subMatrixName){
    cout << "\t";
    for(int i = 0; i < subMatrixName.length(); i++){
        cout << subMatrixName[i] << "\t";
    }
    cout << endl;
    cout << endl;
    for(int i = 0; i < arr.size(); i++){
        cout << subMatrixName[i] << "\t";
        for(int j = 0; j < arr.size(); j++){
            cout << arr[i][j] << "\t";
        }
        cout << endl;
    }
}

int determineLocation(char letter){
    int index;
    switch (letter){
    case 'A':
        index = 0;
        break;
    case 'R':
        index = 1;
        break;
    case 'N':
        index = 2;
        break;
    case 'D':
        index = 3;
        break;
    case 'C':
        index = 4;
        break;
    case 'Q':
        index = 5;
        break;
    case 'E':
        index = 6;
        break;
    case 'G':
        index = 7;
        break;
    case 'H':
        index = 8;
        break;
    case 'I':
        index = 9;
        break;
    case 'L':
        index = 10;
        break;
    case 'K':
        index = 11;
        break;
    case 'M':
        index = 12;
        break;
    case 'F':
        index = 13;
        break;
    case 'P':
        index = 14;
        break;
    case 'S':
        index = 15;
        break;
    case 'T':
        index = 16;
        break;
    case 'W':
        index = 17;
        break;
    case 'Y':
        index = 18;
        break;
    case 'V':
        index = 19;
        break;
    default:
        break;
    }

    return index;
}

vector<vector<float>> globalAlignment(string gene1, string gene2, vector<vector<float>> subMatrix, float gapPenalty, string& optPath){
    int column = gene1.size();
    int row = gene2.size();
    string path;


    // Initialize the OPT and Path Table
    vector<vector<float>> OPT(row + 1, vector<float>(column + 1, 0));
    vector<vector<char>> Path(row + 1, vector<char>(column + 1));

    for(int i = 0; i <= column; i++){
        OPT[0][i] = i * gapPenalty;
        Path[0][i] = 'Z';
    }

    for(int j = 0; j <= row; j++){
        OPT[j][0] = j * gapPenalty;
        Path[j][0] = 'Z';
    }
    Path[0][0] = 'Z';
    
    OPT[1][1] = max({subMatrix[0][0] + OPT[1-1][1-1], gapPenalty + OPT[1-1][1], gapPenalty + OPT[1][1-1]});

    
    // Fill in the OPT and Path table
     for(int i = 1; i <= row; i++){
        for(int j = 1; j <= column; j++){
            if(gene2[i-1] == gene1[j - 1]){
                OPT[i][j] = max({subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1], gapPenalty + OPT[i - 1][j], gapPenalty + OPT[i][j - 1]});
                if(OPT[i][j] == subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1]) Path[i][j] = 'D';
                else if(OPT[i][j] == gapPenalty + OPT[i - 1][j]) Path[i][j] = 'V';
                else Path[i][j] = 'H';
            }else{
                OPT[i][j] = max({subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1], gapPenalty + OPT[i - 1][j], gapPenalty + OPT[i][j - 1]});
                if(OPT[i][j] == subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1]) Path[i][j] = 'D';
                else if(OPT[i][j] == gapPenalty + OPT[i - 1][j]) Path[i][j] = 'V';
                else Path[i][j] = 'H';
            }
        }
        
    }

    // Traveling through the Path table to find route of optimal path 
    while(Path[row][column] != 'Z'){
        path += Path[row][column];
        if(Path[row][column] == 'D'){
            row = row - 1;
            column = column - 1;
        }else if(Path[row][column] == 'V'){
            row = row - 1;
        }else{
            column = column - 1;
        }
    }

    optPath = path;

    return OPT;
}

vector<vector<float>> localAlignment(string gene1, string gene2, vector<vector<float>> subMatrix, float gapPenalty, string& optPath, string& alignedGene1, string& alignedGene2, int& rowMax, int& columnMax){
    int column = gene1.size(), row = gene2.size(), currentMax = 0;
    string path;

    // Initialize the OPT table
    vector<vector<float>> OPT(row + 1, vector<float>(column + 1, 0));
    vector<vector<char>> Path(row + 1, vector<char>(column + 1));

    for(int i = 0; i <= column; i++){
        OPT[0][i] = 0;
        Path[0][i] = 'Z';
    }

    for(int j = 0; j <= row; j++){
        OPT[j][0] = 0;
        Path[j][0] = 'Z';
    }
    
     for(int i = 1; i <= row; i++){
        for(int j = 1; j <= column; j++){
            if(gene2[i-1] == gene1[j - 1]){
                OPT[i][j] = max({subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1], gapPenalty + OPT[i - 1][j], gapPenalty + OPT[i][j - 1]});
                if(OPT[i][j] >= currentMax){
                    rowMax = i;
                    columnMax = j;
                    currentMax = OPT[i][j];
                }
                if(OPT[i][j] == subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1]) Path[i][j] = 'D';
                else if(OPT[i][j] == gapPenalty + OPT[i - 1][j]) Path[i][j] = 'V';
                else Path[i][j] = 'H';
                if(OPT[i][j] < 0){
                    OPT[i][j] = 0;
                    Path[i][j] = 'Z';
                }
            }else{
                OPT[i][j] = max({subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1], gapPenalty + OPT[i - 1][j], gapPenalty + OPT[i][j - 1]});
                if(OPT[i][j] >= currentMax){
                    rowMax = i;
                    columnMax = j;
                    currentMax = OPT[i][j];
                }
                if(OPT[i][j] == subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1]) Path[i][j] = 'D';
                else if(OPT[i][j] == gapPenalty + OPT[i - 1][j]) Path[i][j] = 'V';
                else Path[i][j] = 'H';
                if(OPT[i][j] < 0){
                    OPT[i][j] = 0;
                    Path[i][j] = 'Z';
                }
            }
        }
        
    }

    row = rowMax;
    column = columnMax;

    while(Path[row][column] != 'Z'){
        path += Path[row][column];
        if(Path[row][column] == 'D'){
            row = row - 1;
            column = column - 1;
        }else if(Path[row][column] == 'V'){
            row = row - 1;
        }else{
            column = column - 1;
        }
    }

    optPath = path;

  // Generate the aligned sequences based on the optimal path
    row = rowMax;
    column = columnMax;
    while(Path[row][column] != 'Z'){
        if(Path[row][column] == 'D'){
            alignedGene1 += gene1[column - 1];
            alignedGene2 += gene2[row - 1];
            row--;
            column--;
        } else if(Path[row][column] == 'V'){
            alignedGene1 += '-';
            alignedGene2 += gene2[row - 1];
            row--;
        } else if(Path[row][column] == 'H'){
            alignedGene1 += gene1[column - 1];
            alignedGene2 += '-';
            column--;
        }
    }

    return OPT;   
}

vector<vector<float>> semiAlignment(string gene1, string gene2, vector<vector<float>> subMatrix, float gapPenalty, string& optPath, string& alignedGene1, string& alignedGene2, int& rowMax, int& columnMax){
    int column = gene1.size(), row = gene2.size(), currentMax = 0;
    string path;

    // Initialize the OPT table
    vector<vector<float>> OPT(row + 1, vector<float>(column + 1, 0));
    vector<vector<char>> Path(row + 1, vector<char>(column + 1));

    for(int i = 0; i <= column; i++){
        OPT[0][i] = 0;
        Path[0][i] = 'Z';
    }

    for(int j = 0; j <= row; j++){
        OPT[j][0] = 0;
        Path[j][0] = 'Z';
    }
    
     for(int i = 1; i <= row; i++){
        for(int j = 1; j <= column; j++){
            if(gene2[i-1] == gene1[j - 1]){
                OPT[i][j] = max({subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1], gapPenalty + OPT[i - 1][j], gapPenalty + OPT[i][j - 1]});
                if(OPT[i][j] >= currentMax){
                    rowMax = i;
                    columnMax = j;
                    currentMax = OPT[i][j];
                }
                if(OPT[i][j] == subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1]) Path[i][j] = 'D';
                else if(OPT[i][j] == gapPenalty + OPT[i - 1][j]) Path[i][j] = 'V';
                else Path[i][j] = 'H';
            }else{
                OPT[i][j] = max({subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1], gapPenalty + OPT[i - 1][j], gapPenalty + OPT[i][j - 1]});
                if(OPT[i][j] >= currentMax){
                    rowMax = i;
                    columnMax = j;
                    currentMax = OPT[i][j];
                }
                if(OPT[i][j] == subMatrix[determineLocation(gene2[i - 1])][determineLocation(gene1[j -1])] + OPT[i -1][j - 1]) Path[i][j] = 'D';
                else if(OPT[i][j] == gapPenalty + OPT[i - 1][j]) Path[i][j] = 'V';
                else Path[i][j] = 'H';
            }
        }
        
    }

    row = rowMax;
    column = columnMax;

    while(Path[row][column] != 'Z'){
        path += Path[row][column];
        if(Path[row][column] == 'D'){
            row = row - 1;
            column = column - 1;
        }else if(Path[row][column] == 'V'){
            row = row - 1;
        }else{
            column = column - 1;
        }
    }

    optPath = path;

  // Generate the aligned sequences based on the optimal path
    row = rowMax;
    column = columnMax;
    while(Path[row][column] != 'Z'){
        if(Path[row][column] == 'D'){
            alignedGene1 += gene1[column - 1];
            alignedGene2 += gene2[row - 1];
            row--;
            column--;
        } else if(Path[row][column] == 'V'){
            alignedGene1 += '-';
            alignedGene2 += gene2[row - 1];
            row--;
        } else if(Path[row][column] == 'H'){
            alignedGene1 += gene1[column - 1];
            alignedGene2 += '-';
            column--;
        }
    }

    reverse(alignedGene1.begin(), alignedGene1.end());
    reverse(alignedGene2.begin(), alignedGene2.end());
   
    // Add parentheses around the parts of the words not used in the alignment
    alignedGene1 = "(" + gene1.substr(0, column) + ")" + alignedGene1 + "(" + gene1.substr(columnMax) + ")";
    alignedGene2 = "(" + gene2.substr(0, row) + ")" + alignedGene2 + "(" + gene2.substr(rowMax) + ")";

    return OPT;   
}

void printAlignedSequence(string gene1, string gene2, string path, string type){
    if(type == "global"){
        reverse(path.begin(), path.end());

        if(gene1.length() < gene2.length()){
            while(path.length() < gene2.length()){
                path += 'Z';
            }
        }else{
            while(path.length() < gene1.length()){
                path += 'Z';
            }
        }

        reverse(path.begin(), path.end());

        if(gene2.length() < gene1.length()){
            for(int i = 0; i < gene1.length(); i++){
                cout << gene1[i] << " ";
            }
            cout << endl;
            int counter = 0;
            for(int i = 0; i < path.length(); i++){
                if((path[i] == 'D' || path[i] == 'X') && counter <= path.length()){
                    cout << gene2[counter] << " ";
                    counter++;
                }else{
                    cout << "- ";
                }
            }
            cout << endl;
        }else{
            for(int i = 0; i < gene2.length(); i++){
                cout << gene2[i] << " ";
            }
            cout << endl;
            int counter = 0;
            for(int i = 0; i < gene2.length(); i++){
                if((path[i] == 'D' || path[i] == 'X') && counter <= path.length()){
                    cout << gene1[counter] << " ";
                    counter++;
                }else{
                    cout << "- ";
                }
            }
            cout << endl;
        }
    }else if (type == "local"){
        reverse(gene1.begin(), gene1.end());
        reverse(gene2.begin(), gene2.end());

        if(gene1.length() < gene2.length()){
            for(int i = 0; i < gene2.length(); i++){
                cout << gene2[i] << " ";
            }
            cout << endl;

            for(int i = 0; i < gene1.length(); i++){
                cout << gene1[i] << " ";
            }
        }else{
            for(int i = 0; i < gene1.length(); i++){
                cout << gene1[i] << " ";
            }
            cout << endl;
            for(int i = 0; i < gene2.length(); i++){
                cout << gene2[i] << " ";
            }
        }
    }else{
        cout << gene1 << endl << gene2 << endl;
    }

    cout << endl;

}

