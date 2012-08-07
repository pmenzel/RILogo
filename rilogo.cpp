/*************************************************

  RILogo

  Author: Peter Menzel <pmenzel@gmail.com> and
          Stefan Seemann <seemann@rth.dk>

  Copyright 2012 Peter Menzel

 	This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program, see file LICENCE.
  If not, see <http://www.gnu.org/licenses/>.
  
  See the file README for documentation.

**************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>
#include <map>
#include <math.h>
#include <getopt.h>
#include <unistd.h>
#include <algorithm>
#include <vector>
#include <exception>

#include "ConfigFile.h"

#define MAX(a,b) ((a)>(b)?(a):(b))

#define MIN3(x,y,z)  ((y) <= (z) ? \
                         ((x) <= (y) ? (x) : (y)) \
                     : \
                         ((x) <= (z) ? (x) : (z)))

#define MAX3(x,y,z)  ((y) >= (z) ? \
                         ((x) >= (y) ? (x) : (y)) \
                     : \
                         ((x) >= (z) ? (x) : (z)))


struct rgb_color {
    float r, g, b;    /* Channel intensities between 0.0 and 1.0 */
};

struct hsv_color {
    float h;        /* Hue degree between 0.0 and 360.0 */
    float s;        /* Saturation between 0.0 (gray) and 1.0 */
    float v;        /* Value between 0.0 (black) and 1.0 */
};


struct Config {
	rgb_color arcs_color_from;
	rgb_color arcs_color_to;
	rgb_color arcs_color_nomi;
	rgb_color color_A;
	rgb_color color_C;
	rgb_color color_G;
	rgb_color color_U;
	rgb_color color_M;
	bool debug;
	bool debugSVG;
	bool log_colors;
	bool no_M;
	bool draw_legend_mi;
};


struct rgb_color hex2rgb(std::string hex_string);
std::string rgb2css(struct rgb_color);
struct hsv_color hsv_gradient(hsv_color from, hsv_color to, float percentage);
struct rgb_color hsv_to_rgb(struct hsv_color hsv);
struct hsv_color rgb_to_hsv(struct rgb_color rgb);
void findPairs(const char *ss, const int len, std::map<int,int> *basepairs, char open, char close);
bool findPairsIA(const char *ss, const int len, int middle, std::map<int,int> *basepairs, char open, char close);
int findHighestArc(std::map<int,int> &bp);
void readinput(std::istream & is, std::string & seq, std::map<std::string,std::string> * seqs, std::string & ia, int * profile[], int & numseqs, bool reverse);
void readinputonefile(std::istream & is, std::string & seq1, std::string & seq2, std::map<std::string,std::string> * seqs1, std::map<std::string,std::string> * seqs2, std::string & ia1, std::string & ia2, int * profile1[], int *profile2[], int & numseqs);
std::string trim(std::string& s,const std::string& drop);
int char2index(char c);
char index2char(int i);
std::string getColor(Config * config, float i); // calculate color, where i is a percentage (0.0 - 1.0) on the gradient
void drawLogo(Config * config, int ** profile, float * mi, int length, int numseqs, bool norm_treedist);
void drawLogoScale(Config * config, bool norm_treedist);
void drawRuler(int length, bool reverse, int x_offset, int y_offset,bool no_dots);
void usage(char *progname);
template <class T> void freeProfile(T ** profile);

void mutual_IC(int ** profile, std::map<int,int> *basepairs, std::map<std::string,std::string> * seqs, int length, int numseqs, float * mi, bool reverse, string mi_name);
void mutual_IC_1to2(int ** profile1, int ** profile2, std::map<int,int> *basepairs, std::map<std::string,std::string> * seqs1, std::map<std::string,std::string> * seqs2, int length1, int length2, int numseqs, float * mi, string mi_name);
void treebased_MIC(float ** profile, std::map<int,int> *basepairs, std::map<std::string,std::string> * seqs, int length, int numseqs, float * mi, std::map<std::string,float> * weight, float treedistsum, bool reverse, string mi_name);
void treebased_MIC_1to2(float ** profile1, float ** profile2, std::map<int,int> *basepairs, std::map<std::string,std::string> * seqs1, std::map<std::string,std::string> * seqs2, int length1, int length2, int numseqs, float * miia, std::map<std::string,float> * weight, float treedistsum, string mi_name);

void unionorganism(std::map<std::string,std::string> * seqs1, std::map<std::string,std::string> * seqs2, std::vector<std::string> & unionorg);
void getprofile(std::map<std::string,std::string> * seqs, int length, std::vector<string> & org, std::map<std::string,float> * treedist, float * profile[], bool reverse);
void initweight(std::map<std::string,std::string> * seqs1, std::map<std::string,std::string> * seqs2, std::map<std::string,float> * weight);
float readtreedistfile(std::istream & treedistfile, std::map<std::string,float> * treedist, std::vector<std::string> & unionorg);


using namespace std;

int main(int argc, char **argv) {
	string seq1 = "";
	string seq2 = "";
	string ia1 =  "";
	string ia2 =  "";
	map<string,string> seqs1;
	map<string,string> seqs2;
	int **profile1= new int*[5];
	int **profile2= new int*[5];
	int numseqs1 = 0;  //number of seqs in first alignment
	int numseqs2 = 0;  //         ...     second alignment
	bool draw_logos1 = false;
	bool draw_logos2 = false;
	int len_seq1;
	int len_seq2;
	vector<string> unionorg;
	float **unionorgprofile1 = new float*[5];
	float **unionorgprofile2 = new float*[5];
	bool inter_IA = false;  // will be set to true, if intermolecular base pairs are in the structure
	string filename1 = "";
	string filename2 = "";
	string configfilename = "";
	bool twoinputfiles = false;
	string treedistfilename = "";
	bool use_tree = false;
	float treedistsum = 0;
	//float max_mi = 2.0; // MI, MI^WP, treeMI, treeMI^WP are all max 2bit
	//float min_mi = 0.0;

	// Variables for drawing the SVG
	int margin_text_v = 1; // vertical margin above and below text
	int text_font_size = 12; 
	int text_space1 = text_font_size-3; 
	int text_space2 = text_font_size-3; 
	int legend = 18;
	int arc_pushdown = 2;
	int ruler_space = 6;
	const string logo_font = "Verdana";
	int space_right = 20;
	int space_left = 35;

	// Modifiable options
	bool no_chars = false; // when true, then there will be no logos or ACGT letter printed 
	string mi_name = "miwp";  // either mi or miwp
	bool onlycheckinput = false;
	bool norm_treedist = false;
	bool draw_rulers = true;
	bool verbose = false; // when true, print the MI values to stderr

  Config * config = new Config;
	config->debug = false;
	config->debugSVG = false;
	config->log_colors = true;
	config->no_M = false;
	config->draw_legend_mi = true;
	
	// hard coded defaults, can be overriden by config file
	config->arcs_color_from = hex2rgb("#00ff00");
	config->arcs_color_to = hex2rgb("#ff0000");
	config->color_A = hex2rgb("#14e01a");
	config->color_C = hex2rgb("#007fff");
	config->color_G = hex2rgb("#ff6321");
	config->color_U = hex2rgb("#b7140b");
	config->color_M = hex2rgb("#07570b");
	config->arcs_color_nomi = hex2rgb("#0261c8");


	// --------------------- START ------------------------------------------------------------------
	
	// Read command line params
	char c;
	while ((c = getopt (argc, argv, "hdgzwvc:t:m:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'd':
				config->debug = true; break;
			case 'g':
				config->debugSVG = true; break;
			case 'c':
				configfilename = optarg; break;
			case 't':
				use_tree = true; treedistfilename = optarg; break;
			case 'm':
				mi_name = optarg; break;
			case 'z':
				onlycheckinput = true; break;
			case 'w':
				norm_treedist = true; break;
			case 'v':
				verbose = true; break;
			default:
				usage(argv[0]);
		}
	}

	if(!(mi_name == "mi" || mi_name == "miwp")) { usage(argv[0]); }

	if(!use_tree) norm_treedist = true; // to make sure to use normal 2.0 scaling when no treeMI is used

	if(optind < argc) { // is there an additional parameter for file name?
		filename1 = argv[optind];
		if(optind +1 < argc) {  // there is a second input file name parameter
			filename2 = argv[optind+1];
			twoinputfiles = true;
		}
	}

	try {
		if(configfilename.length() > 0) {
			ConfigFile cf(configfilename,"=","%","" );
			string col;
			cf.readInto<string>(col, "arcs_color_from"); 
			if(col.length()>0)  config->arcs_color_from = hex2rgb(col);
			col = "";
			cf.readInto<string>(col, "arcs_color_to"); 
			if(col.length()>0)  config->arcs_color_to = hex2rgb(col);
			col = "";
			cf.readInto<string>(col, "color_A"); 
			if(col.length()>0)  config->color_A = hex2rgb(col);
			col = "";
			cf.readInto<string>(col, "color_C");
			if(col.length()>0)  config->color_C = hex2rgb(col);
			col = "";
			cf.readInto<string>(col, "color_G");
			if(col.length()>0)  config->color_G = hex2rgb(col);
			col = "";
			cf.readInto<string>(col, "color_U");
			if(col.length()>0)  config->color_U = hex2rgb(col);
			col = "";
			cf.readInto<string>(col, "color_M");
			if(col.length()>0)  config->color_M = hex2rgb(col);

			cf.readInto<bool>(draw_rulers, "draw_rulers");
			cf.readInto<bool>(no_chars, "no_chars");
			cf.readInto<bool>(config->log_colors, "log_colors");
			cf.readInto<bool>(config->no_M, "no_M");
			cf.readInto<bool>(config->draw_legend_mi, "draw_legend_mi");

		}
	}
	catch( ConfigFile::file_not_found& e ) {
		cerr << "Error - File '" << e.filename << "' not found." << endl;
	}

	// read input file
	if(filename1.length() > 0) {
		ifstream inputfile1;
		inputfile1.open(filename1.c_str());
		if(!inputfile1.is_open()) {
 	    cerr << "Could not open file " << filename1 << endl;
			exit(EXIT_FAILURE);
 	  }
		if(!twoinputfiles) {
			readinputonefile(inputfile1,seq1,seq2,&seqs1,&seqs2,ia1,ia2,profile1,profile2,numseqs1);
			numseqs2 = numseqs1;
			inputfile1.close();
			for(map<string,string>::iterator it=seqs1.begin(); it!=seqs1.end(); ++it)
				unionorg.push_back(it->first);
		}
		else {
			ifstream inputfile2;
			inputfile2.open(filename2.c_str());
			if(!inputfile2.is_open()) {
				cerr << "Could not open file " << filename2 << endl;
				exit(EXIT_FAILURE);
			}
			readinput(inputfile1,seq1,&seqs1,ia1,profile1,numseqs1,false);
			inputfile1.close();
			readinput(inputfile2,seq2,&seqs2,ia2,profile2,numseqs2,true);
			inputfile2.close();

			// get profiles for union of organisms in both alignments
			unionorganism(&seqs1, &seqs2, unionorg);
			
		}
	}
	else {  // try to read from STDIN
		readinputonefile(cin,seq1,seq2,&seqs1,&seqs2,ia1,ia2,profile1,profile2,numseqs1);
		numseqs2 = numseqs1;
		for(map<string,string>::iterator it=seqs1.begin(); it!=seqs1.end(); ++it)
			unionorg.push_back(it->first);
	}

  len_seq1 = ia1.length();
  len_seq2 = ia2.length();

	if(config->debug) {
		for(vector<string>::iterator it=unionorg.begin();it!=unionorg.end();it++)
			cerr << " " << *it;
		cerr << endl;
	}

	// check if multiple sequences are there, then switch to logo mode
	if(numseqs1 > 1) { 
		draw_logos1 = true;
		text_space1 = 30;
	}
	if(numseqs2 > 1) { 
		draw_logos2 = true;
		text_space2 = 30;
	}

	map<int, int> bp1;  // store base pairs  1->12, 2->11, etc.. also reverse 12->1,...
	map<int, int> bp2;
	map<int, int> bp1to2; //intermolecular bp.

	// 1. parse IA strings and make interaction map for both intramolecular IAs
	// find basic intramolecular IAs
	findPairs(ia1.c_str(), len_seq1, &bp1, '<', '>');
	findPairs(ia1.c_str(), len_seq1, &bp1, '(', ')');
	findPairs(ia1.c_str(), len_seq1, &bp1, '{', '}');
	findPairs(ia1.c_str(), len_seq1, &bp1, '[', ']');
	findPairs(ia2.c_str(), len_seq2, &bp2, '<', '>');
	findPairs(ia2.c_str(), len_seq2, &bp2, '(', ')');
	findPairs(ia2.c_str(), len_seq2, &bp2, '{', '}');
	findPairs(ia2.c_str(), len_seq2, &bp2, '[', ']');

	// find pseudoknotted pairs, using rfam annotation, uppercase char is 5' end of bp..
	for(char c='a'; c<='z'; c++) {
		char upper = toupper(c);
		// if c or upper(c) be found in both seqs, then it's intermolecular IA, otherwise just normal bp
		if(string::npos != ia1.find_first_of(upper)) { // found uppercase char in 1st ia
			if(string::npos != ia2.find_first_of(c)) { // found corresponding lowercase char in 2nd ia
				// check that lowercase is not in 1st and uppercase is not in 2nd
				if(string::npos != ia1.find_first_of(c)) {  
					cerr << "Found character " << c << " in the structure annotation for both sequences." << endl;
					freeProfile<int>(profile1); freeProfile<int>(profile2);
					exit(EXIT_FAILURE);
				}
				if(string::npos != ia2.find_first_of(upper)) { 
					cerr << "Found character " << upper << " in the structure annotation for both sequences." << endl;
					freeProfile<int>(profile1); freeProfile<int>(profile2);
					exit(EXIT_FAILURE);
				}

				//  map from 1st to 2nd sequence:
				inter_IA = inter_IA | findPairsIA((ia1+ia2).c_str(), len_seq1+len_seq2, len_seq1 , &bp1to2, upper, c);
			}
			else { // only in first
				findPairs(ia1.c_str(), len_seq1, &bp1, upper, c);
			}
		}
		else if(string::npos != ia2.find_first_of(c) || string::npos != ia2.find_first_of(upper)) { // only in 2nd
			findPairs(ia2.c_str(), len_seq2, &bp2, upper, c);
		}
		// else char not found anywhere
	}


	// store weight, e.g. average distance in a tree, of each organism
	// that will be used for calculating the mutual information content
	map<string,float> weight;
	// initialize weight for the organisms in joint of both alignments with 1
	initweight(&seqs1, &seqs2, &weight);

	float * mi1 = NULL, * mi2 = NULL, * miia = NULL;
	bool use_mi1 = false, use_mi2 = false, use_miia = false;  //flags to indicate that MI is calculated, i.e. that there is a MSA as input


	// read average distance of each leaf (organism) in a tree to all other organisms
	if(use_tree) {
		ifstream treedistfile;
		treedistfile.open(treedistfilename.c_str());
		if(!treedistfile.is_open()) {
			cerr << "Could not open file " << treedistfile << endl;
			exit(EXIT_FAILURE);
		}

		try {
			treedistsum = readtreedistfile(treedistfile, &weight, unionorg);
		}
		catch(string & s ) {
			cerr << "The sequence with name " << s << " has no entry in the tree file." << endl;
			treedistfile.close();
			freeProfile<int>(profile1); freeProfile<int>(profile2); freeProfile<float>(unionorgprofile1); freeProfile<float>(unionorgprofile2);
			exit(EXIT_FAILURE);
		}
		treedistfile.close();

		if(!norm_treedist) treedistsum=1.;
		/*
		if(config->debug)
			for(map<string,float>::iterator it=weight.begin(); it!=weight.end(); ++it)
				cerr << it->first << "\t" << it->second << endl;
		*/

		getprofile(&seqs1, ia1.length(), unionorg, &weight, unionorgprofile1, false);
		getprofile(&seqs2, ia2.length(), unionorg, &weight, unionorgprofile2, true);
		/*if(config->debug)
			for(int i=0; i<5; i++) {
				for(unsigned int j=0; j<ia1.length(); j++)
					cerr << " " << unionorgprofile1[i][j];
				cerr << "  ";
				for(unsigned int j=0; j<ia2.length(); j++)
					cerr << " " << unionorgprofile2[i][j];
				cerr << endl;
			}
			*/
	}


	// all input is processed, now can start calculating MI and print the SVG
	// or exit if only the input needed to be checked
	if(onlycheckinput) { 
		freeProfile<int>(profile1);
		freeProfile<int>(profile2);
		freeProfile<float>(unionorgprofile1);
		freeProfile<float>(unionorgprofile2);
		return EXIT_SUCCESS;
	}

	/* --------- calculate the Mutual Information ---------------------------- */
	if(numseqs1 > 1) {
		use_mi1 = true;
		mi1 = new float[len_seq1];
		if(use_tree) {
			treebased_MIC(unionorgprofile1, &bp1, &seqs1, len_seq1, unionorg.size(), mi1, &weight, treedistsum, false, mi_name);
		}
		else {
			mutual_IC(profile1, &bp1, &seqs1, len_seq1, numseqs1, mi1, false, mi_name);
		}
		if(verbose) { cerr << "MI for 1st alignment:"<< endl; for(int i=0; i<len_seq1; i++) if(mi1[i] > 0.0) fprintf(stderr, "col=%i MI=%1.4f\n",i+1, mi1[i]*2); }
	}

	if(numseqs2 > 1) {
		use_mi2 = true;
		mi2 = new float[len_seq2]; // mi2 gets filled reversed
		if(use_tree) {
			treebased_MIC(unionorgprofile2, &bp2, &seqs2, len_seq2, unionorg.size(), mi2, &weight, treedistsum, true, mi_name);
		}
		else {
			mutual_IC(profile2, &bp2, &seqs2, len_seq2, numseqs2, mi2, true, mi_name);
		}
		if(verbose) { cerr << "MI for 2nd alignment:"<< endl; for(int i=0; i<len_seq2; i++) if(mi2[len_seq2 - i -1] > 0.0) fprintf(stderr, "col=%i MI=%1.4f\n",i+1, mi2[len_seq2 - i -1]*2); }
	}

	if(numseqs1 > 1 && numseqs2 > 1) {
		miia = new float[len_seq1];
		use_miia = true;
		if(use_tree) {
			treebased_MIC_1to2(unionorgprofile1, unionorgprofile2, &bp1to2, &seqs1, &seqs2, len_seq1, len_seq2, unionorg.size(), miia, &weight, treedistsum, mi_name);
		}
		else {
			mutual_IC_1to2(profile1, profile2, &bp1to2, &seqs1, &seqs2, len_seq1, len_seq2, numseqs1, miia, mi_name);
		}
		if(verbose) { cerr << "MI for interaction:"<< endl; for(int i=0; i<len_seq1; i++)  if(miia[i] > 0.0)  fprintf(stderr, "col1=%i col2=%i MI=%1.4f\n",i+1, bp1to2[i]+1, miia[i]*2);  }

		// add the MI from the interaction to the other two arrays 
		// remember that miwp2 is reversed
		for(int i=0; i < len_seq1; i++) {
			if(miia[i] > 0.0) {
				mi1[i] = miia[i];
				mi2[len_seq2 - bp1to2[i] - 1] = miia[i];
			}
		}
	}

	if(!use_mi1 && !use_mi2) {
		config->arcs_color_from = config->arcs_color_nomi;
	}
	
	int offset1=0, offset2=0;	
	// 3. calculate the "optimal" x_offset for the 2nd sequence, if offset is <0,
	// then the 1st sequence will be shifted to the right by this amount
	// The optimum is determined by minimizing the sum of squared x-distances between
	// the two bases in an intermoleculare base pair.
	if(inter_IA) {
		int lsm_offset = 0;
		int lsm = INT_MAX;
		for(int curr_offset = -1 * max(len_seq2,len_seq1); curr_offset < max(len_seq2,len_seq1); curr_offset++) {
			//calc square mean for curr_offset
			int sm = 0;
			for(map<int,int>::iterator it	= bp1to2.begin(); it != bp1to2.end(); it++) {
				if(bp1to2.count((*it).first) == 1 && bp1to2[(*it).first] >= 0) {
					int pos1 = (*it).first;
					int pos2 = bp1to2[(*it).first];
					int dist = abs(pos1 - (curr_offset + len_seq2 - pos2) +1); 
					//cerr << "dist from " << pos1 << " to " << pos2 << " (" << curr_offset + len_seq2 - pos2 << ") is " << dist << endl;
					sm += dist * dist;
				}
			}
			if(sm < lsm) {
				lsm_offset = curr_offset;
				lsm = sm;
			}
		}

		if(lsm_offset >= 0) {
			offset1 = 0;
			offset2 = lsm_offset * 10;
		}
		else{
			offset1 = -1 *lsm_offset *10;
			offset2 = 0;
		}
	}
	
	// 4. Draw something Goddamnit!!
	unsigned int distance_between_seq = (inter_IA) ? 40 : 20;
	unsigned int margin = 5;

	int highest_arc1 = (findHighestArc(bp1) * 5 + 5) / arc_pushdown + margin; // find the max of each set of arcs
	int highest_arc2 = (findHighestArc(bp2) * 5 + 5) / arc_pushdown + margin;

	unsigned int y1t = highest_arc1 + margin_text_v ; //top of the first seq
	unsigned int y1c = y1t + margin_text_v + text_space1;  // bottom of text
	unsigned int y1b = y1c + margin_text_v + ruler_space+ margin_text_v; // where intermolecular IA start
	unsigned int y2t = y1b + distance_between_seq;
	unsigned int y2r = y2t + ruler_space ; // where ruler starts
	unsigned int y2c = y2r + margin_text_v + text_space2;  // bottom of text
	unsigned int y2b = y2c + margin_text_v;

	if(no_chars) {
		y1b = y1c = y1t;
		y2b = y2c = y2t;
		draw_logos1 = false;
		draw_logos2 = false;
	} 

	unsigned int height= y2b + highest_arc2 + legend;
	unsigned int width = space_left + space_right + max(len_seq1,len_seq2) * 10  +max(abs(offset1),abs(offset2)); 
	cout << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height << "\">" << endl;
  cout << "<defs>" << endl;
	cout << "<style type=\"text/css\"><![CDATA[" << endl;
	cout << "  .outline{ fill:none; stroke:#0077bb; stroke-width:1 }" << endl;
	cout << "  .fontstyle{ font-family:Verdana; font-size:"<<text_font_size<<"px }" << endl;
	cout << "  .fontstyle53{ fill:#666; font-family:Verdana; font-size:9px }" << endl;
	cout << "  .fontstylescale{ fill:#666; font-family:Verdana; font-size:5px }" << endl;
	cout << "  .fontstyleruler{ fill:#666; font-family:Verdana; font-size:6px; text-anchor:middle;}" << endl;
	cout << "  .fontstylelogo{ font-family:"<<logo_font<<"; font-size:"<<text_font_size<<"px }" << endl;
	cout << "  .logoA{ fill:"<< rgb2css(config->color_A) << ";}" << endl;
	cout << "  .logoC{ fill:"<< rgb2css(config->color_C) << ";}" << endl;
	cout << "  .logoG{ fill:"<< rgb2css(config->color_G) << ";}" << endl;
	cout << "  .logoU{ fill:"<< rgb2css(config->color_U) << ";}" << endl;
	cout << "  .logoM{ fill:"<< rgb2css(config->color_M) << ";}" << endl;
	cout << "]]>" << endl;
  cout << "</style>" << endl;
	if(config->draw_legend_mi && (use_mi1 || use_mi2)) {
		cout << "<linearGradient id=\"grad1\" x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">" << endl;
		cout <<	"<stop offset=\"0%\" style=\"stop-color:"<< getColor(config, 0.0) <<";stop-opacity:1\" />" << endl;
		cout << "<stop offset=\"50%\" style=\"stop-color:"<< getColor(config, 0.5) <<";stop-opacity:1\" />" << endl; //NOTE: only one color in the middle!!
		cout << "<stop offset=\"100%\" style=\"stop-color:"<< getColor(config, 1.0) << ";stop-opacity:1\" />" << endl;
		cout <<	"</linearGradient>" << endl;
	}
  cout << "</defs>" << endl;
  cout << "<rect x=\"0\" y=\"0\" width=\"" << width <<"\" height=\"" << height <<"\" style=\"fill:white;\"/>" << endl;

	cout << "<g transform=\"translate(" << space_left << ")\">" << endl;

	// 1st sequence
	// 5'---------------------------------3'
	// draw intramolecular IAs 
	if(config->debugSVG) {
		cout <<	"<line x1=\""<< offset1 << "\" y1=\""<<y1t<<"\" x2=\""<< offset1 + len_seq1*10  <<"\" y2=\""<<y1t<<"\" stroke=\"#ff0000\"/>" << endl;  
		cout <<	"<line x1=\""<< offset1 << "\" y1=\""<<y1c<<"\" x2=\""<< offset1 + len_seq1*10  <<"\" y2=\""<<y1c<<"\" stroke=\"#ff0000\"/>" << endl;  
		cout <<	"<line x1=\""<< offset1 << "\" y1=\""<<y1b<<"\" x2=\""<< offset1 + len_seq1*10  <<"\" y2=\""<<y1b<<"\" stroke=\"#ff0000\"/>" << endl;  
	}
	if(draw_logos1) {
		cout << "<g transform=\"translate("<< offset1 << " " << y1t  <<")\">" << endl;
			cout << "<g transform=\"translate("<< -2 << " " << 0  <<")\">" << endl;
				cout << "<text class=\"fontstyle53\" y=\"18\" x=\"-28\">5'</text>" << endl;
				cout << "<text class=\"fontstyle53\" y=\"18\" x=\""<< len_seq1*10 + 5 <<"\">3'</text>" << endl;
				drawLogoScale(config,norm_treedist);
			cout << "</g>" << endl;
			drawLogo(config,profile1, mi1,ia1.length(),numseqs1, norm_treedist);
		cout << "</g>" << endl;
	}
	else { // 5' and 3' labels
			int y = (no_chars) ? 5 : 9;
			cout << "<text class=\"fontstyle53\" y=\""<< y1t + y <<"\" x=\""<< offset1 -12 <<"\">5'</text>" << endl;
			cout << "<text class=\"fontstyle53\" y=\""<< y1t + y <<"\" x=\""<< offset1 + len_seq1*10 + 3 <<"\">3'</text>" << endl;
	}
	for(int i = 0; i < len_seq1; i++) {
		if(!draw_logos1) {
			if(no_chars) {  // draw circles
				cout << "<circle cx=\"" <<  offset1 + i * 10 +4  << "\" cy=\"" << y1c << "\" r=\"1\" stroke=\"black\" stroke-width=\"1\" fill=\"black\"/>" << endl;
			}
			else { // draw letters
				cout << "<g transform=\"translate("<< offset1 + i * 10 << " " << y1c  <<")\">" << endl;
				cout << "<text class=\"fontstyle logo"<<seq1[i] <<"\" transform=\"scale(1 1)\">" << seq1[i] << "</text>" << endl;
				cout << "</g>" << endl;
			}
		}

		if(bp1.count(i)==1 && bp1[i] >= 0 && i < bp1[i] ) { // pos i forms bp and i denotes the 5'-base of the pair, then draw an arc
				unsigned int left_x = i * 10 + 5;
				unsigned int right_x = bp1[i] * 10 + 5;
				unsigned int dist = right_x - left_x;
				//unsigned int center = left_x + dist / 2;
				unsigned int radius = dist / 2;
				// MI goes up to 2.0, mi1[] stores half of the MI value, so goes up to 1.0, which can be used as the percentage for  getColor()
				cout << "<path style=\"stroke:" << getColor(config,  use_mi1 ? mi1[i] : 0)  << "\" class=\"outline\" d=\"M " << offset1 + left_x << "," << highest_arc1 << " A" <<radius << "," << radius/arc_pushdown << " 0 " << 1 << "," << 1 << " " << offset1 + right_x << "," << highest_arc1 <<"\"/>" << endl;
		}
	}
	// second sequence
	// 3'---------------------------------5'
	int right_side = len_seq2*10-10;
	if(config->debugSVG) {
		cout <<	"<line x1=\""<< offset2 <<"\" y1=\""<<y2t<<"\" x2=\""<< offset2 + right_side + 10 <<"\" y2=\""<<y2t<<"\" stroke=\"#ff0000\"/>" << endl;  
		cout <<	"<line x1=\""<< offset2 <<"\" y1=\""<<y2r<<"\" x2=\""<< offset2 + right_side + 10 <<"\" y2=\""<<y2r<<"\" stroke=\"#ff0000\"/>" << endl;  
		cout <<	"<line x1=\""<< offset2 <<"\" y1=\""<<y2c<<"\" x2=\""<< offset2 + right_side + 10 <<"\" y2=\""<<y2c<<"\" stroke=\"#ff0000\"/>" << endl;  
		cout <<	"<line x1=\""<< offset2 <<"\" y1=\""<<y2b<<"\" x2=\""<< offset2 + right_side + 10 <<"\" y2=\""<<y2b<<"\" stroke=\"#ff0000\"/>" << endl;  
	}
	if(draw_logos2) {
		cout << "<g transform=\"translate("<< offset2 << " " << y2r  <<")\">" << endl;
			cout << "<g transform=\"translate("<< -2 << " " << 0  <<")\">" << endl;
				cout << "<text class=\"fontstyle53\" y=\"18\" x=\"-28\">3'</text>" << endl;
				cout << "<text class=\"fontstyle53\" y=\"18\" x=\""<< right_side + 15 <<"\">5'</text>" << endl;
				drawLogoScale(config,norm_treedist);
			cout << "</g>" << endl;
			drawLogo(config,profile2, mi2,ia2.length(),numseqs2, norm_treedist);
		cout << "</g>" << endl;
	}
	else { // 5' and 3' labels
			int y = (no_chars) ? 0 : 8;
			cout << "<text class=\"fontstyle53\" y=\""<< y2r + y <<"\" x=\""<< offset2 -12 <<"\">3'</text>" << endl;
			cout << "<text class=\"fontstyle53\" y=\""<< y2r + y <<"\" x=\""<< offset2 + len_seq2*10 + 3 <<"\">5'</text>" << endl;
	}
	// draw intramolecular IAs 
	for(int i = 0; i < len_seq2; i++) {
		 // draw letter
		if(!draw_logos2) {
			if(no_chars) { // draw circles
				cout << "<circle cx=\"" <<  offset2 + right_side - i * 10 +4  << "\" cy=\"" << y2c << "\" r=\"1\" stroke=\"black\" stroke-width=\"1\" fill=\"black\"/>" << endl;
			} else {
				cout << "<text class=\"fontstyle logo"<< seq2[i] <<"\" x=\""<< offset2 + right_side - i * 10 << "\" y=\"" << y2c << "\">" << seq2[i] << "</text>" << endl;
			}
		}
		if(bp2.count(i)==1 && bp2[i] >= 0 && i < bp2[i] ) { // pos i forms bp and i denotes the 5'-base of the pair, then draw an arc
				//convert x coordinates so that they would start from the right side
				unsigned int left_x = right_side - (bp2[i] -1) * 10 - 5;
				unsigned int right_x = right_side - i * 10 + 5 ;
				unsigned int dist = right_x - left_x;
				//unsigned int center = left_x + dist / 2;
				unsigned int radius = dist / 2;
				cout << "<path style=\"stroke:"<< getColor(config, use_mi2 ? mi2[len_seq2 - i -1] : 0) <<"\" class=\"outline\" d=\"M " << offset2 + left_x << "," << y2b << " A" <<radius << "," << radius/arc_pushdown << " 0 " << 0 << "," << 0 << " " << offset2 + right_x << "," << y2b <<"\"/>" << endl;
		}
	}

	if(draw_rulers) {
		if(no_chars) {
			drawRuler(ia1.length(),false,offset1,y1b+6,true);
			drawRuler(ia2.length(),true,offset2,y2r,true);
		}
		else {
			drawRuler(ia1.length(),false,offset1,y1b-2,false);
			drawRuler(ia2.length(),true,offset2,y2r-2,false);
		}
	}

	// draw intermolecular IAs
	if(inter_IA) {
		for(int i = 0; i < len_seq1; i++) {
			if(bp1to2.count(i)==1 && bp1to2[i] >= 0) { // pos i forms IA
					unsigned int left_x = i * 10 + 5 -1;
					unsigned int right_x = right_side - bp1to2[i] * 10 + 5 -1;
					cout <<	"<line style=\"stroke:"<< getColor(config, use_miia ? miia[i] : 0) << "\" x1=\""<< offset1 + left_x << "\" y1=\""<<y1b<<"\" x2=\""<< offset2 + right_x  <<"\" y2=\""<<y2t<<"\" stroke=\"#ff0000\"/>" << endl;  
			}
		}
	}

	if(config->draw_legend_mi && (use_mi1 || use_mi2)) {
		// draw a color gradient legend for mutual information
		cout << "<text class=\"fontstyle\" x=\""<< -12  << "\" y=\"" << height - 9 << "\">" << "MI" << "</text>" << endl;
		cout << "<rect x=\"" << 20   <<"\" y=\"" << height - legend << "\" width=\"60\" height=\"10\" rx=\"1\" ry=\"1\" style=\"fill:url(#grad1);stroke:#000000;stroke-width:1;\"/>" << endl;
		cout << "<text class=\"fontstyleruler\" style=\"text-anchor:start;\" x=\""<< 20  << "\" y=\"" << height - 2  << "\">" << "0" << "</text>" << endl;
		cout << "<text class=\"fontstyleruler\" style=\"text-anchor:end;\" x=\""<< 80  << "\" y=\"" << height - 2 << "\">" << "2.0" << "</text>" << endl;
		cout << "<text class=\"fontstyleruler\" style=\"text-anchor:start;\" x=\""<< -10  << "\" y=\"" << height - 2 << "\">" << "bits" << "</text>" << endl;
	}

	cout << "</g>" << endl;
	cout << "</svg>" << endl;	


	freeProfile<int>(profile1); freeProfile<int>(profile2);
	if(use_mi1)
		free(mi1);
	if(use_mi2)
		free(mi2);
	if(use_miia)
		free(miia);

	return EXIT_SUCCESS;

}

void drawRuler(int length, bool reverse, int x_offset, int y_offset, bool no_dots) {
	int y = y_offset;
	int width = 10;
	if(!reverse)
		for(int i = 0; i < length; i++) {
			if(i==0 || (i+1)%10==0 || i+1 == length) { //number
				cout << "<text class=\"fontstyleruler\" x=\""<< x_offset +  width*i + 4  << "\" y=\"" << y + 1 <<"\">"<< i+1 <<"</text>" << endl;
			}
			else if(!no_dots) {
				if((i+1)%5==0) { //medium line
					cout <<	"<line x1=\""<< x_offset + width*i + 4 << "\" y1=\""<< y+1 <<"\" x2=\""<< x_offset + width*i +4  <<"\" y2=\""<< y-2 <<"\" stroke=\"#444\"/>" << endl;  
				}
				else { //tiny line
					cout <<	"<line x1=\""<< x_offset + width*i + 4 << "\" y1=\""<< y+1 <<"\" x2=\""<< x_offset + width*i +4  <<"\" y2=\""<< y <<"\" stroke=\"#444\"/>" << endl;  
				}
			}
		}
	else
		for(int i = length-1; i >= 0; i--) {
			if(i==length-1 || (length - i)%10==0 || i == 0) { //number
				cout << "<text class=\"fontstyleruler\" x=\""<<  x_offset + width*i + 4  << "\" y=\"" << y + 1 << "\">"<< length - i <<"</text>" << endl;
			}
			else if(!no_dots) {
				if((length - i)%5==0) { //medium line
					cout <<	"<line x1=\""<<x_offset +  width*i + 4 << "\" y1=\""<< y <<"\" x2=\""<<x_offset +  width*i +4  <<"\" y2=\""<< y-3 <<"\" stroke=\"#444\"/>" << endl;  
				}
				else { //tiny line
					cout <<	"<line x1=\""<<x_offset +  width*i + 4 << "\" y1=\""<< y-2 <<"\" x2=\""<<x_offset +  width*i +4  <<"\" y2=\""<< y-3 <<"\" stroke=\"#444\"/>" << endl;  
				}
			}
		}
}

void drawLogo(Config * config, int ** profile, float * mi, int length, int numseqs, bool norm_treedist) {
	float column_height = 2.0;
	if(!config->no_M && !norm_treedist) column_height =  2.2; // set this to 2.2 if normal weighting is used, since the letter stack + M >= 2bit
	float bit_max = 2.0;
	//int num_possible_chars = 4;
	float x = 0;
	int counts [4];
	int indizes [4];

	for(int n = 0; n < length; n++) {
		float y = 30;
		//calculate overall height of the stack
		float bit_height = 0.0; 
		for(int i = 0; i < 4; i++) {
			float freq = (float)profile[i][n] / (float)numseqs;
			bit_height += (freq == 0.0) ? 0.0 : freq * log(freq) / log(bit_max);
			//cerr << bit_height  << endl;
			counts[i] = profile[i][n];
			indizes[i] = i;
		}
		bit_height *= -1;
		float bit_diff = bit_max - bit_height; 
		//cerr << bit_diff << endl; //bit_diff is overall height of stack

		//simple bubble sort of order of bases on the stack, smallest base goes to bottom
		bool sorted = false;
		while(!sorted) {
			sorted = true;
			for(int i = 0; i < 3; i++) {
				if(counts[i] > counts[i+1]) {
					int buf = counts[i+1];
					counts[i+1] = counts[i];
					counts[i] = buf;
					buf = indizes[i+1];
					indizes[i+1] = indizes[i];
					indizes[i] = buf;
					sorted = false;
				}
			}
		}

		for(int i = 0; i < 4; i++) {
			float freq = counts[i]  / (float)numseqs;
			float char_height_px = freq * bit_diff * ( 30 / column_height);
      float scale_factor = char_height_px / 8.8;

			if(scale_factor > 0.0) {
				cout << "<g transform=\"translate("<< x << " " <<  + y  <<")\">" << endl;
				cout << "<text x=\"0\" y=\"0\" class=\"fontstylelogo logo"<< index2char(indizes[i]) <<"\" transform=\"scale(1,"<< scale_factor <<")\">" << index2char(indizes[i]) << "</text>" << endl;
				cout << "</g>" << endl;
			}

			y -= char_height_px;
		}

		// draw an M on top
		if(!config->no_M && mi[n] > 0.0) {
			//cerr << mi[n] << endl;	
			
			float char_height_px = mi[n] * ( 30 / column_height);
      float scale_factor = char_height_px / 8.8;

			if(scale_factor > 0.0) {
				cout << "<g transform=\"translate("<< x << " " <<  + y  <<")\">" << endl;
				cout << "<text x=\"0\" y=\"0\" class=\"fontstylelogo logoM\" transform=\"scale(1,"<< scale_factor <<")\">" << "M" << "</text>" << endl;
				cout << "</g>" << endl;
			}

		}



		x += 10;
	}
}

void drawLogoScale(Config * config, bool norm_treedist) {
		string toplabel = "2.0";
		if(!norm_treedist && !config->no_M) toplabel = "2.2";
		string middlelabel = norm_treedist ? "1.0" : "1.1";
		cout << "<line width=\"1px\" stroke=\"#666\" x1=\"-3\" x2=\"-3\" y1=\""<< 0 << "\" y2=\""<< 30  <<"\"/>" << endl;
		cout << "<line width=\"1px\" stroke=\"#666\" x1=\"-5\" x2=\"-3\" y1=\""<< 0 << "\" y2=\""<< 0  <<"\"/>" << endl;
		cout << "<line width=\"1px\" stroke=\"#666\" x1=\"-5\" x2=\"-3\" y1=\""<< 15 << "\" y2=\""<< 15  <<"\"/>" << endl;
		cout << "<line width=\"1px\" stroke=\"#666\" x1=\"-5\" x2=\"-3\" y1=\""<< 30 << "\" y2=\""<< 30  <<"\"/>" << endl;
		cout << "<text class=\"fontstylescale\" x=\"-15\" y=\"32\">0.0</text>" << endl;
		cout << "<text class=\"fontstylescale\" x=\"-15\" y=\"2\">"<< toplabel <<"</text>" << endl;
		cout << "<text class=\"fontstylescale\" x=\"-15\" y=\"17\">"<< middlelabel <<"</text>" << endl;
		cout << "<text class=\"fontstylescale\" x=\"-15\" y=\"-5\">bits</text>" << endl;
}

int findHighestArc(std::map<int,int> &bp) {

	map<int,int>::iterator it;
	int highest_diff = 0;
	for(it=bp.begin(); it != bp.end(); it++) {
		if( (*it).second >= 0) {
			highest_diff = max(highest_diff,  (abs((*it).second) - abs((*it).first)));
		}	
	}
	return highest_diff;
}

void readinputonefile(istream & is, std::string & seq1, std::string & seq2, std::map<std::string,std::string> * seqs1, std::map<std::string,std::string> * seqs2, std::string & ia1, std::string & ia2, int * profile1[], int *profile2[],int & numseqs) {

	// profile2 gets filled reversed!!

	// read fasta
	string line;
	string ia_line = "";
	string seq_line = "";
	bool read_seq = false;
	bool read_structure = false;
	string curr_name = "";
	std::map<std::string,std::string> seqs;

	while(!is.eof()) {
		getline(is, line);
		if(line[0]=='>') {  // start new entry
			read_structure = false;
			read_seq = false;
			if(seq_line.length() > 0) {  // store old sequence
				seqs.insert(make_pair(curr_name,seq_line));
				seq_line = "";
				curr_name = "";
			}

			if(line.find("structure") != string::npos) {  //this is the structure annotation line
				read_structure = true;	
			}
			else { //normal sequence line
				read_seq = true;
				string n = line.substr(1);
				curr_name = trim(n, " \r\t\n\v\f\b");
			}
		}
		else if(read_structure) {
			//append to structure line
			ia_line += trim(line, " \r\t\n\v\f\b");
		}
		else if(read_seq) {
			seq_line += trim(line," \r\t\n\v\f\b");
		}
	}

	numseqs = seqs.size();
	if(numseqs == 0) {
		cerr << "No sequences were found." << endl;
		exit(EXIT_FAILURE);
	}

	if(ia_line.length() == 0) {
		cerr << "No structure annotation line was found." << endl;
		exit(EXIT_FAILURE);
	}
	// parse ia line
	// should only contain [A-Za-z&.()<>-]+
	size_t found = ia_line.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVQXYZabcdefghijklmnopqrstuvwxyz&()<>{}[]-.");
	if(string::npos != found) {
		cerr << "Found invalid character " << ia_line[found] << " at position " << found << " in the structure line!" << endl;
		exit(EXIT_FAILURE);
	}
	//check if exactly one & is present in this line
	found = ia_line.find_first_of("&");
	if(string::npos == found) {  // not found
		cerr << "Separator & was not found in the structure line!" << endl;
		exit(EXIT_FAILURE);
	}
	else {  // found &, look for another one
		if(string::npos != ia_line.find_first_of("&",found+1)) {
			cerr << "Found more than one separators & in the structure line" << endl;
			exit(EXIT_FAILURE);
		}
	}
	// split at & 
	ia1 = ia_line.substr(0,found);
	ia2 = ia_line.substr(found+1);

	for(int i = 0 ; i < 5 ; i++)
		profile1[i] = new int[ia1.length()];
	for(int i = 0 ; i < 5 ; i++)
		for(unsigned int j = 0 ; j < ia1.length(); j++)
			profile1[i][j] = 0;

	for(int i = 0 ; i < 5 ; i++ )
		profile2[i] = new int[ia2.length()];
	for(int i = 0 ; i < 5 ; i++ )
		for(unsigned int j = 0 ; j < ia2.length(); j++)
			profile2[i][j] = 0;

	// parse each sequence line in the seqs map 
	for(map<string,string>::iterator it = seqs.begin(); it != seqs.end(); ++it) {
		// should only contain [A-Za-z&]+, no gaps anyway
		string seq_line = it->second;
		size_t found = seq_line.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVQXYZabcdefghijklmnopqrstuvwxyz-.&");
		if(string::npos != found) {
			cerr << "Found invalid character " << seq_line[found] << " at position " << found << " in the sequence "<< it->first << endl;
			freeProfile<int>(profile1); freeProfile<int>(profile2);
			exit(EXIT_FAILURE);
		}
		//check if exactly one & is present in this line
		found = seq_line.find_first_of("&");
		if(string::npos == found) {  // not found
			cerr << "Separator & was not found in the sequence "<< it->first << endl;
			freeProfile<int>(profile1); freeProfile<int>(profile2);
			exit(EXIT_FAILURE);
		}
		else {  // found &, look for another one
			if(string::npos != seq_line.find_first_of("&",found+1)) {
				cerr << "Found more than one separator & in the sequence line" << endl;
				freeProfile<int>(profile1); freeProfile<int>(profile2);
				exit(EXIT_FAILURE);
			}
		}
		// split at & 
		seq1 = seq_line.substr(0,found);
		seq2 = seq_line.substr(found+1);

		//check for same lengths
		if(seq1.length() != ia1.length()) {
			cerr << "The length of the first sequence ("<< seq1.length()  <<") for " << it->first << " and the length of the structure annotation (" <<  ia1.length() << ") are not the same." << endl;
			freeProfile<int>(profile1); freeProfile<int>(profile2);
			exit(EXIT_FAILURE);
		}
		if(seq2.length() != ia2.length()) {
			cerr << "The length of the second sequence ("<< seq2.length()  <<") for " << it->first <<" and the length of the structure annotation (" <<  ia2.length() << ") are not the same." << endl;
			freeProfile<int>(profile1); freeProfile<int>(profile2);
			exit(EXIT_FAILURE);
		}
	
		//add seq1 to profile1
		for(unsigned int i = 0; i < ia1.length(); i++) {
			char b = seq1[i];
			profile1[char2index(b)][i]++;
		}
		//add seq2 to profile2, reversed
		for(unsigned int i = 0; i < ia2.length(); i++) {
			char b = seq2[ia2.length() - i -1];
			profile2[char2index(b)][i]++;
		}

		//add seq1 to seqs1
		seqs1->insert(make_pair(it->first,seq1));
		//add seq2 to seqs2
		seqs2->insert(make_pair(it->first,seq2));
	} // end for each sequence 

}

void readinput(std::istream & is, std::string & seq, std::map<std::string,std::string> * seqs, std::string & ia, int * profile[], int & numseqs, bool reverse) {

	// read fasta
	string line;
	string ia_line = "";
	string seq_line = "";
	bool read_seq = false;
	bool read_structure = false;
	string curr_name = "";

	while(!is.eof()) {
		getline(is, line);
		if(line[0]=='>') {  // start new entry
			read_structure = false;
			read_seq = false;
			if(seq_line.length() > 0) {  // store old sequence
				seqs->insert(make_pair(curr_name,seq_line));
				seq_line = "";
				curr_name = "";
			}

			if(line.find("structure") != string::npos) {  //this is the structure annotation line
				read_structure = true;	
			}
			else { //normal sequence line
				read_seq = true;
				string n = line.substr(1);
				curr_name = trim(n, " \r\t\n\v\f\b");
			}
		}
		else if(read_structure) {
			//append to structure line
			ia_line += trim(line, " \r\t\n\v\f\b");
		}
		else if(read_seq) {
			seq_line += trim(line," \r\t\n\v\f\b");
		}
	}

	numseqs = seqs->size();
	if(numseqs == 0) {
		cerr << "No sequences were found." << endl;
		exit(EXIT_FAILURE);
	}

	if(ia_line.length() == 0) {
		cerr << "No structure annotation line was found." << endl;
		exit(EXIT_FAILURE);
	}
	// parse ia line
	// should only contain [A-Za-z.()<>-]+
	size_t found = ia_line.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVQXYZabcdefghijklmnopqrstuvwxyz()<>{}[]-.");
	if(string::npos != found) {
		cerr << "Found invalid character " << ia_line[found] << " at position " << found << " in the structure line!" << endl;
		exit(EXIT_FAILURE);
	}
	ia = ia_line;

	for(int i = 0 ; i < 5 ; i++)
		profile[i] = new int[ia.length()];
	for(int i = 0 ; i < 5 ; i++)
		for(unsigned int j = 0 ; j < ia.length(); j++)
			profile[i][j] = 0;

	// parse each sequence line in the seqs map 
	//for(map<string,string>::iterator it = seqs->begin(); it != seqs->end(); ++it) {
	for(map<string,string>::iterator it = seqs->begin(); it != seqs->end(); ++it) {
		string seq_line = it->second;
		size_t found = seq_line.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVQXYZabcdefghijklmnopqrstuvwxyz-.");
		if(string::npos != found) {
			cerr << "Found invalid character " << seq_line[found] << " at position " << found << " in the sequence "<< it->first << endl;
			freeProfile<int>(profile);
			exit(EXIT_FAILURE);
		}
		seq = seq_line;

		//check for same lengths
		if(seq.length() != ia.length()) {
			cerr << "The length of the sequence ("<< seq.length()  <<") for " << it->first << " and the length of the structure annotation (" <<  ia.length() << ") are not the same." << endl;
			freeProfile<int>(profile);
			exit(EXIT_FAILURE);
		}
	
		if(!reverse) //add seq to profile in normal direction
			for(unsigned int i = 0; i < ia.length(); i++) {
				char b = seq[i];
				profile[char2index(b)][i]++;
			}
		else //add seq to profile, reversed
			for(unsigned int i = 0; i < ia.length(); i++) {
				char b = seq[ia.length() - i -1];
				profile[char2index(b)][i]++;
			}

	} // end for each sequence 


}


void findPairs(const char *ss, const int len, map<int,int> *basepairs, char open, char close) {

	//cerr << "finding pairs using " << open << " and " << close << " in seq " << ss << endl;
	const char single = '.';
	const char gap = '-';
	stack<int> lager;

	for(int i = 0; i < len; i++) {

		// This is the left side of a base pair
		if(ss[i] == open) {
			lager.push(i);
		}
		else if(ss[i] == close) { // right side of a bp

			if(lager.empty()) {
				cerr << "ERROR: Secondary structure annotation is not well-formed. (Underflow at position " << i+1 <<")"<< endl;
				exit(EXIT_FAILURE);
			}
			int start = lager.top(); // get top element
			lager.pop(); // remove element.. stupid stack impl.
			
			basepairs->insert(make_pair(i, start));
			basepairs->insert(make_pair(start,i));
		}
		else {
			if(ss[i] == single) {
				basepairs->insert(make_pair(i,-1));	
			}
			else if(ss[i] == gap) {
				basepairs->insert(make_pair(i,-1));	
			}
			else {
				//cerr << "Found character >>" << ss[i] << "<< at position " << i+1 << " in the secondary structure string. That might cause problems.. Let's see what happens." << endl;
			}
		}
	}
	
	if(!lager.empty()) {
		cerr << "ERROR: Secondary structure annotation is not well-formed. (Overflow)" << endl;
		exit(EXIT_FAILURE);
	}

}

bool findPairsIA(const char *ss, const int len, int middle, map<int,int> *basepairs, char open, char close) {

	//cerr << "finding pairs using " << open << " and " << close << " in seq " << ss << endl;
	const char single = '.';
	const char gap = '-';
	stack<int> lager;
	bool found_bp = false;
	for(int i = 0; i < len; i++) {

		// This is the left side of a base pair
		if(ss[i] == open) {
			lager.push(i);
		}
		else if(ss[i] == close) { // right side of a bp

			if(lager.empty()) {
				cerr << "ERROR: Secondary structure annotation is not well-formed. (Underflow at position " << i+1 <<")"<< endl;
				exit(EXIT_FAILURE);
			}
			int start = lager.top(); // get top element
			lager.pop(); // remove element.. stupid stack impl.
			
			//basepairs->insert(make_pair(i, start));
			basepairs->insert(make_pair(start,i-middle));
			found_bp = true;
		}
		else {
			if(ss[i] == single) {
				basepairs->insert(make_pair(i,-1));	
			}
			else if(ss[i] == gap) {
				basepairs->insert(make_pair(i,-1));	
			}
			else {
				//LOG(lWARN) << "Found character >>" << ss[i] << "<< at position " << i+1 << " in the secondary structure string. That might cause problems.. Let's see what happens.";
			}
		}
	}
	
	if(!lager.empty()) {
		cerr << "ERROR: Secondary structure annotation is not well-formed. (Overflow)" << endl;
		exit(EXIT_FAILURE);
	}

	return found_bp;

}

void mutual_IC(int ** profile, map<int,int> *basepairs, map<string,string> * seqs, int length, int numseqs, float * mi, bool reverse, string mi_name) {

	bool bgap = true; // this correponds to treeMI^WP
	bool bjan = false;

	if(mi_name ==  "mi") {
		bgap = false;  // == treeMI
		bjan = true;
	}

	// define canonical base pairs: A-U, C-G, G-U, U-A, G-C, U-G
	int canono[6] = {0, 1, 2, 3, 2, 3};
	int canonc[6] = {3, 2, 3, 0, 1, 2};
	// remember which bases (A, C, G, U) occur in a column if base pair
	
	for(int i = 0; i < length; i++)
		mi[i]=0.;

	for(map<int,int>::iterator itcons = basepairs->begin(); itcons != basepairs->end(); itcons++) { // pos i forms bp and i denotes the 5'-base of the pair
		if(itcons->second < 0 || itcons->first > itcons->second)
			continue;

		// mutual information using only canonical basepairs
		// [Lindgreen 2006 -> MI 9: Basepairs + gap penalty]
		char open, close;
		float pBP, relbaseo, relbasec;
		float qBP=0;
		int pair=0;
		bool bpbool;
		float curr_mi = 0.;

		for(int i=0; i<6; i++) {
			bpbool = false;
			for(map<string,string>::iterator j = seqs->begin(); j != seqs->end(); j++) {
				open = (j->second)[itcons->first];
				close = (j->second)[itcons->second];
				//open = (*seqs)[j->first][itcons->first];
				//close = (*seqs)[j->first][itcons->second];
				if( char2index(open) == canono[i] && char2index(close) == canonc[i] ) {
					pair++;
					bpbool = true;
				}
			}
			if(bpbool) {
				if(!reverse) {
					relbaseo = profile[canono[i]][itcons->first];
					relbasec = profile[canonc[i]][itcons->second];
				}
				else {
					relbaseo = profile[canono[i]][length - itcons->first-1];
					relbasec = profile[canonc[i]][length - itcons->second-1];
				}

				qBP += relbaseo*relbasec;
				//fprintf(stderr, "%f %f %f\n", relbaseo,relbasec,qBP);
			}
		}
		pBP = (float)pair/(float)numseqs;
		qBP = qBP/(float)(numseqs*numseqs);
		if(qBP != 0 && pBP != 0)
			curr_mi = pBP*(log(pBP/qBP)/log(2));

		//fprintf(stderr, "%s %i %6.4f %6.4f %6.4f\n",mi_name.c_str(), pair, pBP, qBP, curr_mi);

		// gap penalty [Lindgreen 2006]
		if( bgap ) {
			int gaps = 0;
			for(map<string,string>::iterator j = seqs->begin(); j != seqs->end(); j++) {
				open = (j->second)[itcons->first];
				close = (j->second)[itcons->second];
				if( open=='-' || open=='_' || open=='.' || close=='-' || close=='_' || close=='.' )
					gaps++;
			}
			curr_mi -= (float)gaps / (float)numseqs;
		}
		else // calculate Kullback-Leibler distance between two distributions
			if( bjan )
				if(pBP != 1)
					curr_mi += (1-pBP) * (log( (1-pBP)/(1-qBP) )/log(2));

		if(curr_mi < 0.0) curr_mi = 0.0; 

		if(reverse) {
			mi[length - itcons->first -1] =  curr_mi / 2.0;
			mi[length - itcons->second -1] = mi[length - itcons->first -1];
		}
		else {
			mi[itcons->first] =  curr_mi / 2.0;
			mi[itcons->second] = mi[itcons->first];
		}
	}

	// BEGIN just for statistics
	/*fprintf(stderr,"MI_BP_PenTrueBPs<-c(");
	for(map<int,int>::iterator itcons = basepairs->begin(); itcons != basepairs->end(); itcons++){
		if(itcons->second < 0 || itcons->first > itcons->second)
			continue;
		fprintf(stderr, "%10.9f,", mi[itcons->first]);
	}
	fprintf(stderr,")\n");
	// MIC of non-basepairs
	float nonbpmiwp;
	bool bp;
	fprintf(stderr,"MI_BP_PenNonBPs<-c(");
	for(int i=0;i<length;i++){
		for(int j=i+4;j<length;j++){
			bp = false;
			for(map<int,int>::iterator itcons = basepairs->begin(); itcons != basepairs->end(); itcons++) {
				if(itcons->first==i && itcons->second==j)
					bp = true;
			}
			if(bp) continue;

			char open, close;
			float pBP, relbaseo, relbasec;
			float qBP=0;
			float pair=0;
			bool bpbool;
			map<string,float>::iterator ito;
			for(int k=0; k<6; k++) {
				bpbool = false;
				for(map<string,string>::iterator itseq = seqs->begin(); itseq != seqs->end(); itseq++) {
					open = (itseq->second)[i];
					close = (itseq->second)[j];
					if( char2index(open) == canono[k] && char2index(close) == canonc[k] ) {
						bpbool = true;
						pair++;
					}
				}
				if( bpbool ) {
					if(!reverse) {
						relbaseo = profile[canono[k]][i];
						relbasec = profile[canonc[k]][j];
					}
					else {
						relbaseo = profile[canono[k]][length - i-1];
						relbasec = profile[canonc[k]][length - j-1];
					}
					qBP += relbaseo*relbasec;
				}
			}
			pBP = (float)pair/(float)numseqs;
			qBP = qBP/(float)(numseqs*numseqs);
			if(qBP != 0 && pBP != 0) {
				nonbpmiwp = pBP*(log(pBP/qBP)/log(2));
				int gaps = 0;
				for(map<string,string>::iterator itseq = seqs->begin(); itseq != seqs->end(); itseq++) {
					open = (itseq->second)[i];
					close = (itseq->second)[j];
					if( open=='-' || open=='_' || open=='.' || close=='-' || close=='_' || close=='.' )
						gaps++;
				}
				nonbpmiwp -= (float)gaps / (float)numseqs;
				nonbpmiwp /= 2;
			}
			else {
				nonbpmiwp = 0;
			}
			fprintf(stderr,"%10.9f,",nonbpmiwp);
		}
	}
	fprintf(stderr,")\n");
	// END just for statistics
	 */
}


void treebased_MIC(float ** profile, map<int,int> *basepairs, map<string,string> * seqs, int length, int numseqs, float * mi, map<string,float> * weight, float treedistsum, bool reverse, string mi_name) {

	bool bgap = true; // this correponds to treeMI^WP
	bool bjan = false;

	if(mi_name ==  "mi") {
		bgap = false;  // == treeMI
		bjan = true;
	}

	// define canonical base pairs: A-U, C-G, G-U, U-A, G-C, U-G
	int canono[6] = {0, 1, 2, 3, 2, 3};
	int canonc[6] = {3, 2, 3, 0, 1, 2};

	for(int i = 0; i < length; i++)
		mi[i] = 0.;

	for(map<int,int>::iterator itcons = basepairs->begin(); itcons != basepairs->end(); itcons++) { // pos i forms bp and i denotes the 5'-base of the pair
		if(itcons->second < 0 || itcons->first > itcons->second)
			continue;

		// tree based mutual information
		// using only canonical base pairs (+ gaps) [Lindgreen 2006 -> MI 9]
		char open, close;
		float pBP, relbaseo, relbasec;
		float qBP=0;
		float pair=0;
		bool bpbool;
		string organ;
		float curr_mi = 0.0;
		map<string,float>::iterator ito;

		for(int i=0; i<6; i++) {
			bpbool = false;
			for(map<string,string>::iterator j = seqs->begin(); j != seqs->end(); j++) {
				organ = j->first;
				ito = weight->find(organ);
				if( ito == weight->end() )
					continue;
				open = (j->second)[itcons->first];
				close = (j->second)[itcons->second];
				if( char2index(open) == canono[i] && char2index(close) == canonc[i] ) {
					pair += ito->second;
					bpbool = true;
				}
			}
			if( bpbool ) {
				if(!reverse) {
					relbaseo = profile[canono[i]][itcons->first];
					relbasec = profile[canonc[i]][itcons->second];
				}
				else {
					relbaseo = profile[canono[i]][length - itcons->first-1];
					relbasec = profile[canonc[i]][length - itcons->second-1];
				}
				qBP += relbaseo*relbasec;
			}
		}

		if(treedistsum==1.0) {
			pBP = (float)pair/(float)numseqs;
			qBP = qBP/(float)(numseqs*numseqs);
		}
		else {
			pBP = (float)pair/(float)treedistsum;
			qBP = qBP/(float)(treedistsum*treedistsum);
		}
		if(qBP != 0 && pBP != 0)
			curr_mi = pBP*(log(pBP/qBP)/log(2));

		//fprintf(stderr, "tree%s: %6.4f %6.4f %6.4f %6.4f\n",mi_name.c_str(), pair, pBP, qBP, curr_mi);

		// gap penalty [Lindgreen 2006]
		if( bgap ) {
			int gaps = 0;
			for(map<string,string>::iterator j = seqs->begin(); j != seqs->end(); j++) {
				organ = j->first;
				ito = weight->find(organ);
				if(ito == weight->end())
					continue;
				open = (j->second)[itcons->first];
				close = (j->second)[itcons->second];
				if( open=='-' || open=='_' || open=='.' || close=='-' || close=='_' || close=='.' )
					gaps++;
			}
			curr_mi -= (float)gaps / (float)numseqs;
		}

		// calculate Kullback-Leibler distance between two distributions
		else
			if( bjan )
				if(pBP != 1)
					curr_mi += (1-pBP) * (log( (1-pBP)/(1-qBP) )/log(2));

		if(curr_mi < 0.0) curr_mi = 0.0; 

		if(reverse) {
			mi[length - itcons->first -1] =  curr_mi / 2.0;
			mi[length - itcons->second -1] = mi[length - itcons->first -1];
		}
		else {
			mi[itcons->first] =  curr_mi / 2.0;
			mi[itcons->second] = mi[itcons->first];
		}

	}

	// BEGIN just for statistics
	/*
	fprintf(stderr,"tbMICTrueBPs<-c(");
	for(map<int,int>::iterator itcons = basepairs->begin(); itcons != basepairs->end(); itcons++){
		if(itcons->second < 0 || itcons->first > itcons->second)
			continue;
		fprintf(stderr, "%10.9f,", mit[itcons->first]);
	}
	fprintf(stderr,")\n");
	// tree bases MIC of non-basepairs
	float nonbpmit;
	bool bp;
	fprintf(stderr,"tbMICNonBPs<-c(");
	for(int i=0;i<length;i++){
		for(int j=i+4;j<length;j++){
			bp = false;
			for(map<int,int>::iterator itcons = basepairs->begin(); itcons != basepairs->end(); itcons++) {
				if(itcons->first==i && itcons->second==j)
					bp = true;
			}
			if(bp) continue;

			char open, close;
			float pBP, relbaseo, relbasec;
			float qBP=0;
			float pair=0;
			bool bpbool;
			string organ;
			map<string,float>::iterator ito;
			for(int k=0; k<6; k++) {
				bpbool = false;
				for(map<string,string>::iterator itseq = seqs->begin(); itseq != seqs->end(); itseq++) {
					organ = itseq->first;
					ito = weight->find(organ);
					if( ito == weight->end() )
						continue;
					open = (itseq->second)[i];
					close = (itseq->second)[j];
					if( char2index(open) == canono[k] && char2index(close) == canonc[k] ) {
						pair += ito->second;
						bpbool = true;
					}
				}
				if( bpbool ) {
					if(!reverse) {
						relbaseo = profile[canono[k]][i];
						relbasec = profile[canonc[k]][j];
					}
					else {
						relbaseo = profile[canono[k]][length - i-1];
						relbasec = profile[canonc[k]][length - j-1];
					}
					qBP += relbaseo*relbasec;
				}
			}
			pBP = (float)pair/(float)numseqs;
			qBP = qBP/(float)(numseqs*numseqs);
			if(qBP != 0 && pBP != 0) {
				nonbpmit = pBP*(log(pBP/qBP)/log(2));
				if( bgap ) {
					int gaps = 0;
					for(map<string,string>::iterator itseq = seqs->begin(); itseq != seqs->end(); itseq++) {
						organ = itseq->first;
						ito = weight->find(organ);
						if( ito == weight->end() )
							continue;
						open = (itseq->second)[i];
						close = (itseq->second)[j];
						if( open=='-' || open=='_' || open=='.' || close=='-' || close=='_' || close=='.' )
							gaps++;
					}
					nonbpmit -= (float)gaps / (float)numseqs;
				}
				else
					if( bjan )
						if(pBP != 1)
							nonbpmit += (1-pBP) * (log( (1-pBP)/(1-qBP) )/log(2));
				nonbpmit /= 2;
			}
			else {
				nonbpmit = 0;
			}
			fprintf(stderr,"%10.9f,",nonbpmit);
		}
	}
	fprintf(stderr,")\n");
	*/
	// END just for statistics
}

void mutual_IC_1to2(int ** profile1, int ** profile2, std::map<int,int> *basepairs, std::map<string,string> * seqs1, std::map<string,string> * seqs2, int length1, int length2, int numseqs, float * mi, string mi_name) {

	bool bgap = true; // this correponds to treeMI^WP
	bool bjan = false;

	if(mi_name ==  "mi") {
		bgap = false;  // == treeMI
		bjan = true;
	}

	// define canonical base pairs: A-U, C-G, G-U, U-A, G-C, U-G
	int canono[6] = {0, 1, 2, 3, 2, 3};
	int canonc[6] = {3, 2, 3, 0, 1, 2};
	
	for(int i = 0; i < length1; i++)
		mi[i]=0.;

	for(map<int,int>::iterator itcons = basepairs->begin(); itcons != basepairs->end(); itcons++) { // pos i forms bp and i denotes the 5'-base of the pair
		if(itcons->second < 0)
			continue;
		// mutual information using only canonical basepairs
		// [Lindgreen 2006 -> MI 9: Basepairs + gap penalty]
		char open, close;
		float pBP, relbaseo, relbasec;
		float qBP=0;
		int pair=0;
		bool bpbool;

		for(int i=0; i<6; i++) {
			bpbool = false;
			for(map<string,string>::iterator j = seqs1->begin(); j != seqs1->end(); j++) {
				open = (*seqs1)[j->first][itcons->first];
				close = (*seqs2)[j->first][itcons->second];
				if( char2index(open) == canono[i] && char2index(close) == canonc[i] ) {
					pair++;
					bpbool = true;
				}
			}
			if( bpbool ) {
				relbaseo = profile1[canono[i]][itcons->first];
				relbasec = profile2[canonc[i]][length2 - itcons->second-1];
				qBP += relbaseo*relbasec;
			}
		}
		pBP = (float)pair/(float)numseqs;
		qBP = qBP/(float)(numseqs*numseqs);
		if(qBP != 0 && pBP != 0)
			mi[itcons->first]=pBP*(log(pBP/qBP)/log(2));
		//fprintf(stderr, "%i %6.4f %6.4f %6.4f\n", pair, pBP, qBP, mi[itcons->first]);

		// gap penalty [Lindgreen 2006]
		if( bgap ) {
			int gaps = 0;
			for(map<string,string>::iterator j = seqs1->begin(); j != seqs1->end(); j++) {
				open = (*seqs1)[j->first][itcons->first];
				close = (*seqs2)[j->first][itcons->second];
				if( open=='-' || open=='_' || open=='.' || close=='-' || close=='_' || close=='.' )
					gaps++;
			}
			mi[itcons->first] -= (float)gaps / (float)numseqs;
		}
		else // calculate Kullback-Leibler distance between two distributions
			if( bjan )
				if(pBP != 1)
					mi[itcons->first] += (1-pBP) * (log( (1-pBP)/(1-qBP) )/log(2));

		if(mi[itcons->first] < 0.0) mi[itcons->first] = 0.0; 

		mi[itcons->first] /= 2;
	}
}


void treebased_MIC_1to2(float ** profile1, float ** profile2, std::map<int,int> *basepairs, std::map<string,string> * seqs1, std::map<string,string> * seqs2, int length1, int length2, int numseqs, float * mi, std::map<string,float> * weight, float treedistsum, string mi_name) {

	bool bgap = true; // this correponds to treeMI^WP
	bool bjan = false;

	if(mi_name ==  "mi") {
		bgap = false;  // == treeMI
		bjan = true;
	}

	// define canonical base pairs: A-U, C-G, G-U, U-A, G-C, U-G
	int canono[6] = {0, 1, 2, 3, 2, 3};
	int canonc[6] = {3, 2, 3, 0, 1, 2};

	for(int i = 0; i < length1; i++)
		mi[i]=0.;

	for(map<int,int>::iterator itcons = basepairs->begin(); itcons != basepairs->end(); itcons++) { // pos i forms bp and i denotes the 5'-base of the pair
		if(itcons->second < 0)
			continue;

		// tree based mutual information
		// using only canonical base pairs (+ gaps) [Lindgreen 2006 -> MI 9]
		char open, close;
		float pBP, relbaseo, relbasec;
		float qBP=0;
		float pair=0;
		bool bpbool;
		string organ;
		map<string,float>::iterator ito;

		for(int i=0; i<6; i++) {
			bpbool = false;
			for(map<string,string>::iterator j = seqs1->begin(); j != seqs1->end(); j++) {
				organ = j->first;
				ito = weight->find(organ);
				if( ito == weight->end() )
					continue;
				open = (*seqs1)[j->first][itcons->first];
				close = (*seqs2)[j->first][itcons->second];
				if( char2index(open) == canono[i] && char2index(close) == canonc[i] ) {
					pair += ito->second;
					bpbool = true;
				}
			}
			if( bpbool ) {
				relbaseo = profile1[canono[i]][itcons->first];
				relbasec = profile2[canonc[i]][length2 - itcons->second-1];
				qBP += relbaseo*relbasec;
			}
		}

		if(treedistsum==1.0) {
			pBP = (float)pair/(float)numseqs;
			qBP = qBP/(float)(numseqs*numseqs);
		}
		else {
			pBP = (float)pair/(float)treedistsum;
			qBP = qBP/(float)(treedistsum*treedistsum);
		}

		if(qBP != 0 && pBP != 0)
			mi[itcons->first]=pBP*(log(pBP/qBP)/log(2));

		// gap penalty [Lindgreen 2006]
		if(bgap) {
			int gaps = 0;
			for(map<string,string>::iterator j = seqs1->begin(); j != seqs1->end(); j++) {
				organ = j->first;
				ito = weight->find(organ);
				if( ito == weight->end() )
					continue;
				open = (*seqs1)[j->first][itcons->first];
				close = (*seqs2)[j->first][itcons->second];
				if( open=='-' || open=='_' || open=='.' || close=='-' || close=='_' || close=='.' )
					gaps++;
			}
			mi[itcons->first] -= (float)gaps / (float)numseqs;
		}

		// calculate Kullback-Leibler distance between two distributions
		else
			if( bjan )
				if(pBP != 1)
					mi[itcons->first] += (1-pBP) * (log( (1-pBP)/(1-qBP) )/log(2));

		if(mi[itcons->first] < 0.0) mi[itcons->first] = 0.0;
		mi[itcons->first] /= 2;
		//fprintf(stderr, "mi: %i %6.4f %6.4f %6.4f %6.4f %6.4f\n", numseqs, treedistsum, pair, pBP, qBP, mi[itcons->first]);
	}
}


void unionorganism(std::map<std::string,std::string> * seqs1, std::map<std::string,std::string> * seqs2, std::vector<std::string> & unionorg) {
	for(map<string,string>::iterator it1 = seqs1->begin(); it1 != seqs1->end(); ++it1) {
		for(map<string,string>::iterator it2 = seqs2->begin(); it2 != seqs2->end(); ++it2) {
			if( (it1->first).compare(it2->first) == 0 ) {
				unionorg.push_back(it1->first);
				break;
			}
		}
	}
}


void getprofile(std::map<std::string,std::string> * seqs, int length, std::vector<std::string> & org, std::map<std::string,float> * weight, float * profile[], bool reverse) {

	for(int i = 0 ; i < 5 ; i++)
		profile[i] = new float[length];
	for(int i = 0 ; i < 5 ; i++)
		for(int j = 0 ; j < length; j++)
			profile[i][j] = 0.;

	// count the occurrence of each nucleotide in the sublist of desired organisms
	// and weight it
	for(map<string,string>::iterator it = seqs->begin(); it != seqs->end(); ++it) {
		string organ = "";
		string seq = "";
		for(vector<string>::iterator ito = org.begin(); ito != org.end(); ++ito) {
			if( ito->compare(it->first) == 0 ) {
				organ = it->first;
				seq = it->second;
				break;
			}
		}
		if( seq.length() == 0 )
			continue;

		if(!reverse) //add seq to profile in normal direction
			for(int i = 0; i < length; i++) {
				char b = seq[i];
				profile[char2index(b)][i] += (*weight)[organ];
			}
		else //add seq to profile, reversed
			for(int i = 0; i < length; i++) {
				char b = seq[length - i -1];
				profile[char2index(b)][i] += (*weight)[organ];
			}
	}

}


void initweight(std::map<std::string,std::string> * seqs1, std::map<std::string,std::string> * seqs2, std::map<std::string,float> * weight) {

	for(map<string,string>::iterator it1 = seqs1->begin(); it1 != seqs1->end(); ++it1)
		weight->insert(make_pair(it1->first,1.));

	for(map<string,string>::iterator it2 = seqs2->begin(); it2 != seqs2->end(); ++it2)
		for(map<string,float>::iterator it1 = weight->begin(); it1 != weight->end(); ++it1)
			if( (it2->first).compare(it1->first) != 0 ) {
				weight->insert(make_pair(it2->first,1.));
				break;
			}
}


/*
 * organisms in the the unionorg list get assigned a normalized weight from the average tree distances
 * the other organisms get assigned a zero weight
 */
float readtreedistfile(std::istream & is, std::map<std::string,float> * treedist, std::vector<std::string> & unionorg) {

	string line;
	string name;
	float dist;
	float nd=0;
	float treedistsum=0.;

	for(map<string,float>::iterator it = treedist->begin(); it != treedist->end(); ++it)
			it->second = 0;

	while(!is.eof()) {
		getline(is, line);
		stringstream linestream(line);
		linestream >> name >> dist;
		for(vector<string>::iterator itu=unionorg.begin(); itu!=unionorg.end(); ++itu)
			if( itu->compare(name) == 0 ) {
				if(treedist->count(name)>0) {
					(*treedist)[name] = dist;
					nd += dist;
				}
				break;
			}
		/*		for(map<string,float>::iterator it = treedist->begin(); it != treedist->end(); ++it)
					if( (it->first).compare(name) == 0 ) {
						it->second = dist;
						nd += dist;
						break;
					}*/
	}

	// Go through all names in unionorg and check if they have a tree distance now, otherwise exit
	for(vector<string>::iterator itu=unionorg.begin(); itu!=unionorg.end(); itu++) {
			if(treedist->count(*itu) == 0 || (treedist->count(*itu) == 1 && (*treedist)[*itu]  == 0.0)) throw *itu;
	}

	for(map<string,float>::iterator it = treedist->begin(); it != treedist->end(); ++it)
			if( it->second != 0.0 ) {
				it->second = (nd - it->second) / nd;  // this is the same as 1 - d_avg / nd
				treedistsum += it->second;
			}
	return treedistsum;
}


std::string getColor(Config * config, float value) {
	if(config->log_colors && value > 0.0) {
		//(x+ log(x)/7 + 1)/2)
		value = (value + log(value) / 7 + 1)/2;	
	}
	if(value < 0.0) value = 0.0;  // filter negative log results for veryyyy small values

	struct rgb_color from_rgb = {config->arcs_color_from.r,config->arcs_color_from.g,config->arcs_color_from.b};
	struct rgb_color to_rgb = {config->arcs_color_to.r,config->arcs_color_to.g,config->arcs_color_to.b};
	hsv_color from = rgb_to_hsv(from_rgb);
	hsv_color to = rgb_to_hsv(to_rgb);

	hsv_color final = hsv_gradient(from,to,value);
	rgb_color final_rgb = hsv_to_rgb(final);
	stringstream ostr;
	ostr <<  "rgb(" << round(final_rgb.r * 255) << "," << round(final_rgb.g * 255) << "," << round(final_rgb.b* 255) << ")";
	//if(config->debug) cerr << value << " "  << ostr.str()  << endl;
	return ostr.str();

}


std::string rgb2css(struct rgb_color rgb) {
	stringstream ostr;
	ostr <<  "rgb(" << round(rgb.r * 255) << "," << round(rgb.g * 255) << "," << round(rgb.b* 255) << ")";
	return ostr.str();
}

struct hsv_color hsv_gradient(hsv_color from, hsv_color to, float percentage) { //percentage being between 0 and 1
	// adopted from http://stackoverflow.com/questions/2593832/how-to-interpolate-hue-values-in-hsv-colour-space
	hsv_color final;
	final.s = (1 - percentage) * from.s + percentage * to.s;
	final.v = (1 - percentage) * from.v + percentage * to.v;
	from.h /= 360.0;
	to.h /= 360.0;
	float distCCW = (from.h >= to.h) ? from.h - to.h : 1 + from.h - to.h;
  float distCW  = (from.h >= to.h) ? 1 + to.h - from.h : to.h - from.h;
	float h = (distCW <= distCCW) ? from.h + (distCW * percentage) : from.h - (distCCW * percentage);
	if (h < 0) h = 1 + h;
	if (h > 1) h = h - 1;
	final.h = h*360;

	return final;
}



struct hsv_color rgb_to_hsv(struct rgb_color rgb) {
	 // adopted from http://en.literateprograms.org/RGB_to_HSV_color_space_conversion_%28C%29
    hsv_color hsv;
    float rgb_min, rgb_max;
    rgb_min = MIN3(rgb.r, rgb.g, rgb.b);
    rgb_max = MAX3(rgb.r, rgb.g, rgb.b);

    hsv.v = rgb_max;
    if (hsv.v == 0) {
        hsv.h = hsv.s = 0;
        return hsv;
    }

    /* Normalize value to 1 */
    rgb.r /= hsv.v;
    rgb.g /= hsv.v;
    rgb.b /= hsv.v;
    rgb_min = MIN3(rgb.r, rgb.g, rgb.b);
    rgb_max = MAX3(rgb.r, rgb.g, rgb.b);

    hsv.s = rgb_max - rgb_min;
    if (hsv.s == 0) {
        hsv.h = 0;
        return hsv;
    }

    /* Normalize saturation to 1 */
    rgb.r = (rgb.r - rgb_min)/(rgb_max - rgb_min);
    rgb.g = (rgb.g - rgb_min)/(rgb_max - rgb_min);
    rgb.b = (rgb.b - rgb_min)/(rgb_max - rgb_min);
    rgb_min = MIN3(rgb.r, rgb.g, rgb.b);
    rgb_max = MAX3(rgb.r, rgb.g, rgb.b);

    /* Compute hue */
    if (rgb_max == rgb.r) {
        hsv.h = 0.0 + 60.0*(rgb.g - rgb.b);
        if (hsv.h < 0.0) {
            hsv.h += 360.0;
        }
    } else if (rgb_max == rgb.g) {
        hsv.h = 120.0 + 60.0*(rgb.b - rgb.r);
    } else /* rgb_max == rgb.b */ {
        hsv.h = 240.0 + 60.0*(rgb.r - rgb.g);
    }

    return hsv;

}


struct rgb_color hsv_to_rgb(struct hsv_color hsv) {
	hsv.h /= 60; // sector 0 to 5
	int i = (int)floor( hsv.h );
	float f = hsv.h - i; // factorial part of h

	float p = hsv.v * ( 1 - hsv.s );
	float q = hsv.v * ( 1 - hsv.s * f );
	float t = hsv.v * ( 1 - hsv.s * ( 1 - f ) );
  rgb_color rgb;

	switch( i ) {
		case 0:
             rgb.r = hsv.v;  rgb.g = t;   rgb.b = p;
              break;
		case 1:
              rgb.r = q;   rgb.g = hsv.v;  rgb.b = p;
              break;
		case 2:
              rgb.r = p;   rgb.g = hsv.v;   rgb.b = t;
              break;
		case 3:
              rgb.r = p;   rgb.g = q;   rgb.b = hsv.v;
              break;
		case 4:
              rgb.r = t;   rgb.g = p;   rgb.b = hsv.v;
              break;
		default: // case 5:
              rgb.r = hsv.v;   rgb.g = p;   rgb.b = q;
              break;
	}

	return rgb;
}

struct rgb_color hex2rgb(std::string hex_string) {
	rgb_color rgb;
	rgb.r = 0.0; rgb.b = 0.0; rgb.g = 0.0;

	//cut off first char if it is a #
	if(string::npos != hex_string.find_first_of('#')) { // found # in the string
		hex_string.erase(0,1);
	}
	//check if length is 6
	if(hex_string.length() != 6) {
		cerr << "Could not parse color values from the string \"" << hex_string << "\"! Use a format like #003366" << endl;
		exit(1);
	}
	int i;
	stringstream ss_r(hex_string.substr(0,2));
	if((ss_r >> std::hex >> i).fail()) {
		cerr << "Could not parse color values from the string \"" << hex_string << "\"! Use the hex format like #0033ff" << endl;
		exit(1);
	}
	else { rgb.r = i / 255.0;}

	stringstream ss_g(hex_string.substr(2,2));
	if((ss_g >> std::hex >> i).fail()) {
		cerr << "Could not parse color values from the string \"" << hex_string << "\"! Use the hex format like #0033ff" << endl;
		exit(1);
	}
	else {  rgb.g = i / 255.0; }

	stringstream ss_b(hex_string.substr(4,2));
	if((ss_b >> std::hex >> i).fail()) {
		cerr << "Could not parse color values from the string \"" << hex_string << "\"! Use the hex format like #0033ff" << endl;
		exit(1);
	}
	else {  rgb.b = i / 255.0; }
	return rgb;
}


std::string trim(std::string& s,const std::string& drop = " ") {
		std::string r=s.erase(s.find_last_not_of(drop)+1);
		return r.erase(0,r.find_first_not_of(drop));
}


int char2index(char c) {
	switch (c) {
	    case 'A' : case 'a' : return 0;
	    case 'C' : case 'c' : return 1;
	    case 'G' : case 'g' : return 2;
	    case 'T' : case 't' : case 'U' : case 'u' : return 3;
	    default  : return 4;
	}
}

char index2char(int i) {
	switch (i) {
	    case 0 : return 'A';
	    case 1 : return 'C';
	    case 2 : return 'G';
	    case 3 : return 'U';
			case 5 : return '-';
	    default : return 'N';
	}
}

void usage(char *progname) { 
  fprintf(stderr, "Usage:\n   %s [options] inputfile1.fa [inputfile2.fa]\n", progname);
  fprintf(stderr, "Options are:\n");
  fprintf(stderr, "   -m NAME       Name of mutual information measure. Options are 'mi' and 'miwp'. Default is 'miwp'.\n");
  fprintf(stderr, "   -t FILENAME   Switch to treeMI or treeMI^WP (depending on parameter -m) measure and specify the name of the file with the average tree distances.\n");
  fprintf(stderr, "   -w            Use N_d instead of N in the weighting of observed and expected base pair frequencies.\n");
  fprintf(stderr, "   -c FILENAME   Name of config file\n");
  fprintf(stderr, "   -d            Debug mode\n");
  fprintf(stderr, "   -g            Debug mode for SVG output\n");
  fprintf(stderr, "   -v            Verbose mode, prints calculated MI values to STDERR\n");
	exit(EXIT_FAILURE);
}

template <class T>
void freeProfile(T ** profile){
	for(int i = 0 ; i < 5 ; i++) {
		delete[] profile[i];
	}
	delete[] profile;
}

