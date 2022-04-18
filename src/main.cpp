#include <iostream>       
#include <fstream>
#include <string>        
#include <vector>
#include <algorithm>
#include <bitset>        
#include <cstdlib>
#include <ctime>

const int N=2;
const int K=1;
const int R=4;
const int NB_ERRORS = 2; // Number of errors to add to the transmission channel

//#define DEBUG

using namespace std; 

// Class to model a path for GSM_decoding
class Path {
  public:
    bitset<R> reg;
    vector<bitset<K>> message_decoded;
    int total_dist;

    Path() {};
    Path(bitset<R> a, int b) : reg(a), total_dist(b) {}
};


// Calculating Hamming distance between two bitsets
int hammingDistance(bitset<N> bitset1, bitset<N> bitset2) {
  int cpt = 0;

  if (bitset1[0] != bitset2[0]) {
    cpt++;
  }

  if (bitset1[1] != bitset2[1]) {
    cpt++;
  }

  return cpt; // Number of bits that are different from one bitset to the other
}


// Calculating the 2-bit output from a 5-bit state (convolutionnal coding, like in GSM_code)
bitset<N> code(const bitset<R+1> & reg) {
  bitset<N> output;
	int g0, g1;
  bitset<R+1> G0(25);
  bitset<R+1> G1(27);
	g0 = (reg & G0).count() % 2; // Number of bits set after an AND operation (between initial register and G0 polynom) (0 if even, 1 if odd)
	g1 = (reg & G1).count() % 2; // Same as before but with G1 polynom

  output.reset();
	output.set(0, g0);
	output.set(1, g1);

  return output;
}


// Comparing the cumulated distance of paths
bool distanceComparison(Path path1, Path path2) {
  return ( path1.total_dist < path2.total_dist);
}


////////////////////////////////////////////////////////////
//      template<int bits> bitset<bits> randBitset()      //
//                                                        //
//               Generate random bitset                   //
////////////////////////////////////////////////////////////
template<int bits> bitset<bits> randBitset() { 
  bitset<bits> r(rand());
  for(int i = 0; i < bits/16 - 1; i++)
  {
    r <<= 16;
    r |= bitset<bits>(rand()); 
  }
  return r;
} 


////////////////////////////////////////////////////////////
// vector< bitset<N> > GSM_code(vector< bitset<K> > mess) //
//                                                        //
//     Convolutional coding of a message (GSM norm)       //
////////////////////////////////////////////////////////////
vector< bitset<N> > GSM_code(vector< bitset<K> > mess) {
  int i=0, g0, g1;
  vector< bitset<N> > mess_out;

  bitset<N> cod_out; 
  bitset<R+1> G0(25);
  bitset<R+1> G1(27); 
  bitset<R+1> reg; 
  reg.reset();
  
  #ifdef DEBUG
    cout << "-------------------- Debug Informations (Coding) --------------------" << endl << endl;
    cout << "Initial register ( u(i-4)  u(i-3)  u(i-2)  u(i-1)  u(i)  ): " << reg << endl;
    cout << "Polynom G0       ( g0(i-4) g0(i-3) g0(i-2) g0(i-1) g0(i) ): " << G0 << endl;
    cout << "Polynom G1       ( g1(i-4) g1(i-3) g1(i-2) g1(i-1) g1(i) ): " << G1 << endl << endl;
  #endif

  for (vector<bitset<K> >::iterator it = mess.begin(); it != mess.end(); ++it) {
    reg = reg<<1; // Left shifting of the initial register
    reg.set(0,(*it).count()); // Setting the first bit of the initial register if the current bit scanned is set

    g0 = (reg&G0).count()%2;
    g1 = (reg&G1).count()%2;

    cod_out.reset();
    cod_out.set(0,g0);
    cod_out.set(1,g1);

    mess_out.push_back(cod_out); // Pushing 2-bits binary words at a time in the final message (to output)
    
    #ifdef DEBUG
    cout << "Block number: " << ++i << " - In frame: "<< *it << endl; 
    cout << "\t Current status of registers: "<< reg << endl;
    cout << "\t Out : " << cod_out << endl;
    #endif
  }

  #ifdef DEBUG
    cout << "------------------------------------------------------------------" << endl << endl;
  #endif

  return mess_out;
}


/////////////////////////////////////////////////////////////////////////
// vector< bitset<N> >  GSM_transmission(vector< bitset<N> > mess_cod) //
//                                                                     //
//         Simulation of a transmission channel => adding errors       //
/////////////////////////////////////////////////////////////////////////
vector< bitset<N> >  GSM_transmission(vector< bitset<N> > mess_cod) {
  vector< bitset<N> > mess_tra = mess_cod;

  for (int i = 0; i < NB_ERRORS; i++) { // Injecting the maximum number of errors to the transmission, at random indexes
    int random_index = rand() % mess_cod.size();  // Random N-bit word
    int random_place = rand() % N;                // Random bit in word
    mess_tra[random_index][random_place] = !mess_tra[random_index][random_place]; // Inverting a random bit of the N-bit word designated by the random index
  }
  
  return mess_tra;
}


//////////////////////////////////////////////////////////////////
// vector< bitset<K> > GSM_decode(vector< bitset<N> > mess_tra) //
//                                                              //
//     Convolutional decoding of a message (GSM norm)           //
//////////////////////////////////////////////////////////////////
vector< bitset<K> > GSM_decode(vector< bitset<N> > mess_tra) {
  vector< bitset<K> > mess_dec;

  vector<Path> paths;     // Memorized paths
  vector<Path> new_paths; // Potentially temporary paths

  // Initializing the first path of the graph
  Path initial_path(0000,0);
  paths.push_back(initial_path);

  for (vector<bitset<N>>::iterator it_message = mess_tra.begin() ; it_message != mess_tra.end(); it_message++) {  // For each N-bit word in the transmitted message
    new_paths = paths;
    paths.clear();

    for (vector<Path>::iterator it_path = new_paths.begin() ; it_path != new_paths.end(); it_path++) { // For each memorized path
      // Creating two potential path for the current N-bit word
      Path p1 = *it_path;
      Path p2 = *it_path;

      // Calculating Hamming distance for each potential paths
      bitset<R+1> flow1 (p1.reg.to_string() + "0");                     // Flow = current register + entry bit
      bitset<N> code_flow1 = code(flow1);                               // Convolutionnal coding
      int hamming_dist_1 = hammingDistance(code_flow1, (*it_message));  // Hamming distance calculation
      bitset<R+1> flow2(p2.reg.to_string() + "1");
      bitset<N> code_flow2 = code(flow2);
      int hamming_dist_2 = hammingDistance(code_flow2, (*it_message));

      // Updating the cumulated distance for each path
      p1.total_dist += hamming_dist_1;
      p2.total_dist += hamming_dist_2;

      // Updating their proper register
      p1.reg = (p1.reg<<1); // Shifting
      p2.reg = (p2.reg<<1);
      p1.reg[0] = 0; // Updating first bit value (entry bit)
      p2.reg[0] = 1;

      // Updating the decoded message value for the different paths
      p1.message_decoded.push_back(0);
      p2.message_decoded.push_back(1);

      // Pushing the two potential paths into the global list of paths
      paths.push_back(p1);
      paths.push_back(p2);
    }
  }

  // Sorting the paths based on their cumulated hamming distance
  sort(paths.begin(), paths.end(), distanceComparison);
  mess_dec = paths.at(0).message_decoded; // Selecting the path with the lowest Hamming distance
  
  return mess_dec;
}


//////////////////////////////////////////////////////////////////
//                             MAIN                             //
//////////////////////////////////////////////////////////////////
int main() {
  int NbMot = 12;

  vector< bitset<K> > mess;
  vector< bitset<N> > mess_cod;
  vector< bitset<N> > mess_tra;
  vector< bitset<K> > mess_dec;

  // Random initialization message
  srand( (unsigned)time( NULL ) );
  for(int i=0;i<NbMot;++i)
    mess.push_back(randBitset<K>());
  for(int i=0;i<R;++i)
    mess.push_back(bitset<K>(0));

  // Coding of the message => mess_cod
  mess_cod = GSM_code(mess);

  // Simulation of a transmission (errors) => mess_tra
  mess_tra = GSM_transmission(mess_cod);

  // Decoding of the transmitted message => mess_dec
  mess_dec = GSM_decode(mess_tra);

  cout << "Source Message   : ";
  for (vector<bitset<K> >::iterator it = mess.begin() ; it != mess.end(); ++it)
      cout << ' ' << *it;
    cout << '\n';

  cout << "Coded Message    : ";
  for (vector<bitset<N> >::iterator it = mess_cod.begin() ; it != mess_cod.end(); ++it)
      cout << ' ' << *it;
    cout << '\n';

  cout << "Received Message : ";
  for (vector<bitset<N> >::iterator it = mess_tra.begin() ; it != mess_tra.end(); ++it)
      cout << ' ' << *it;
    cout << '\n';

  cout << "Decoded Message  : ";
  for (vector<bitset<K> >::iterator it = mess_dec.begin() ; it != mess_dec.end(); ++it)
      cout << ' ' << *it;
    cout << '\n';
}

