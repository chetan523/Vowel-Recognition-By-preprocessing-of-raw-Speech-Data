#define _USE_MATH_DEFINES
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#define ll long long
#define ld long double
#define FRAME_SIZE 320
#define p 12
using namespace std;

// This function reads the voice sample file framewise
vector<vector<ld>> readFile(string filename){
	vector<vector<ld>> file;

	ifstream ifs(filename);
	if(!ifs.is_open()){
		cout << "Could not open the file!" << endl;
		system("pause");
		exit(0);
	}
	
	// creating frames of size 320 samples and appending in to file vector
	while(!ifs.eof()){
		vector<ld> frame;
		for(int i=0; i<FRAME_SIZE && !ifs.eof(); ++i){
			ld sample;
			ifs >> sample;
			frame.push_back(sample);
		}

		if(frame.size() < FRAME_SIZE) break;
		file.push_back(frame);
	}

	return file;
}

// This function returns the steady frame of the vowel according maximum STE value frame
vector<vector<ld>> getSteadyFrames(vector<vector<ld>> file, int fileSize){
    ld maxSTE = 0;
    int maxSTE_idx = -1;

    for(int i=0; i<fileSize; ++i){
        ld ste = 0;
        for(int m=0; m<FRAME_SIZE; ++m) ste += file[i][m] * file[i][m];
        
        if(ste > maxSTE){
            maxSTE = ste;
            maxSTE_idx = i;
        }
    }

	// Appending 2 frames before and after the maximum STE frame
    vector<vector<ld>> steadyFrames;
    for(int i=maxSTE_idx-2; i<=maxSTE_idx+2; ++i) steadyFrames.push_back(file[i]);

    return steadyFrames;
}

// This function normalizes the recording file so that maximum amplitutde of the file is 10000
void normalizeValues(vector<vector<ld>>& file, int fileSize){
	ld maxValue = 0;
	for(int i=0; i<fileSize; ++i)
		for(int m=0; m<FRAME_SIZE; ++m)
			maxValue = max(maxValue, abs(file[i][m]));
			
	ld scalingFactor = 10000.0/maxValue;

	for(int i=0; i<fileSize; ++i) 
		for(int m=0; m<FRAME_SIZE; ++m)
			file[i][m] *= scalingFactor;
}

// This function first calculates the DC shift on a silence file and then removes the DC shift from the train and test files
void removeDCShift(vector<vector<ld>>& file, int fileSize){
	ld sum = 0;
	int n = 0;
	ifstream ifs("..\\Sound\\Silence.txt");
	while(!ifs.eof()){
		ld sample;
		ifs >> sample;
		sum += sample;
		++n;
	}

	ld DCShift = sum / (ld)n;

	// removing DC shift from the train and test files
	for(int i=0; i<fileSize; ++i) 
		for(int m=0; m<FRAME_SIZE; ++m)
			file[i][m] -= DCShift;
}

// This function applies the hamming window on a frame
void applyHammingWindow(vector<ld>& frame){
	for(int m=0; m<FRAME_SIZE; ++m){
		ld wm = 0.54 - 0.46 * cos(2*M_PI*m /(FRAME_SIZE-1));
		frame[m] *= wm;
	}
}

// This function calculates the auto correlation terms(r) for the frame
vector<ld> getAutoCorrelationTerms(vector<ld> frame){
	vector<ld> r;
	for(int k=0; k<=p; ++k){
		ld rk = 0;
		for(int m=0; m < FRAME_SIZE-k; ++m)
			rk += frame[m] * frame[m+k];

		r.push_back(rk);
		if(k==0 && rk <= 0) break;
	}
	
	return r;
}

// This function calculates the Linear Predictive Coefficients using the Durbin's Algorithm
vector<ld> applyDurbins(vector<ld> r){
	vector<ld> alpha;
	ld E = r[0];

	for(int i=1; i<=p; ++i){
		ld sum = 0;
		for(int j=1;j<i; ++j) sum += alpha[j] * r[i-j];
		ld k = (r[i] - sum) / E;
		
		vector<ld> alpha_new(i+1);
		alpha_new[i] = k;
		for(int j=1; j<i; ++j) 
			alpha_new[j] = alpha[j] - k * alpha[i-j];
		alpha = alpha_new;
		E = (1 - k*k) * E;
	}

	return alpha;
}

// This function preforms the Linear Predictive Ananlysis and returns the Linear Predictive Coefficients
vector<ld> performLPCAnalysis(vector<ld> frame){
	applyHammingWindow(frame);
	vector<ld> r,a;
	r = getAutoCorrelationTerms(frame);
	if(r[0] <= 0) return a;
	a = applyDurbins(r);
	return a;
}

// This function calculates the Cepstral Coefficients by Linear Transformation of the Linear Predicitve Coefficients
// Note: The Linear Predicitive Coefficients are not inverted for calculating the Cepstral Coefficients
vector<ld> getCepstralCoefficients(vector<ld> a){
	// for(int i=1; i<=p; ++i) a[i] = -a[i];
	vector<ld> C(p+1);
	for(int m=1; m<=p; ++m){
		ld sum = 0;
		for(int k=1; k<m; ++k)
			sum += ((ld)k/(ld)m) * C[k] * a[m-k];
		C[m] = a[m] + sum;
	}

	return C;
}

// This function applies the raised sine window on the cepstral coefficients
void applyRaisedSineWindow(vector<ld>& C){
	for(int m=1; m<=p; ++m){
		ld wm = 1 + (p/2) * sin(M_PI*m/p);
		C[m] *= wm;
	}
}

// This function returns the Cepstrsal Coeffiients calculated on a frame
// Note: Raised sine window is not applied on the cepstral coefficients
vector<vector<ld>> getCoefficients(vector<vector<ld>> steadyFrames){
	vector<vector<ld>> C;
	int n = steadyFrames.size();
	for(int i=0; i<n; ++i){
		vector<ld> lpc;
		lpc = performLPCAnalysis(steadyFrames[i]);

		vector<ld> cc;
		cc = getCepstralCoefficients(lpc);
		// applyRaisedSineWindow(cc);
		C.push_back(cc);
	}

	return C;	
}	

// This function creates the codebook from the training recording files
vector<vector<ld>> train(string file_path){
    char vowels[] = {'a', 'e', 'i', 'o', 'u'};
	vector<vector<ld>> codebook;

    for(int v=0; v<5; ++v){
		vector<vector<vector<ld>>> C;
        for(ll i=1; i<=10; ++i){
            string filename = file_path + "204101041_" + vowels[v] + "_" + to_string(i) + ".txt"; 
			vector<vector<ld>> file;
			file = readFile(filename);
			int fileSize = file.size();

			removeDCShift(file, fileSize);
			normalizeValues(file, fileSize);
			vector<vector<ld>> steadyFrames;
			steadyFrames = getSteadyFrames(file, fileSize);
			vector<vector<ld>> c;
			c = getCoefficients(steadyFrames);
			C.push_back(c);
        }

		
		for(int f=0; f<5; ++f){
			vector<ld> C_avg(p);
			for(int i=0; i<10; ++i){
				for(int j=1; j<=p; ++j){
					C_avg[j-1] += C[i][f][j];
				}
			}
			
			for(int i=0; i<p; ++i) C_avg[i] /= 10;
			codebook.push_back(C_avg);
		}
    }

	return codebook;
}

// This function prints the codebook created using the training recording files
void printCodebook(vector<vector<ld>> codebook){
	char vowels[] = {'a', 'e', 'i', 'o', 'u'};
	for(int f=0; f<25; ++f){
		if(f%5 == 0) printf("\nCode vectors for vowel %c are:\n", vowels[f/5]);
		vector<ld> codeVector = codebook[f];
		for(int i=0; i<p; ++i) cout << codeVector[i] << " ";
		cout << endl;
	}
}

// This function writes the codebook values in a file named "codebook.txt" 
void saveCodebook(vector<vector<ld>> codebook){
	string filename =  "..\\codebook.txt" ;
	ofstream ofs(filename);

	for(int i=0; i<25; ++i){
		for(int j=0; j<p; ++j)
			ofs << codebook[i][j];
		ofs << endl;
	}

	ofs.close();
}

// This function predicts the test file vowel by calculating the Tokhura's distance from all the code vectors of the codebook
char predictVowel(vector<vector<ld>> C, vector<vector<ld>> codebook){
	ld w[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
	char vowels[] = {'a', 'e', 'i', 'o', 'u'};
	vector<ld> TokhuraDistance;

	// Calculating the average Tokhura's Distance
	for(int i=0; i<25; i+=5){
		ld averageTokhuraDistance = 0;
		for(int j=0; j<5; ++j)
			for(int k=0; k<p; ++k)
				averageTokhuraDistance += w[k] * (C[j][k+1] - codebook[j+i][k]) * (C[j][k+1] - codebook[j+i][k]);
			
		averageTokhuraDistance /= 5;
		TokhuraDistance.push_back(averageTokhuraDistance);	
	}

	// Calculating the minimum Tokhura's Distance
	int n = TokhuraDistance.size();
	ld min_dist = 1e300;
	int min_idx = -1;
	for(int i=0; i<n; ++i){
		if(TokhuraDistance[i] < min_dist){
			min_dist = TokhuraDistance[i];
			min_idx = i;
		}
	}

	return vowels[min_idx];
}

// This function predicts the vowel in every test recording file
void test(string file_path, vector<vector<ld>> codebook){
    char vowels[] = {'a', 'e', 'i', 'o', 'u'};
	int correct = 0;
	cout << "\nThe predicted vowels are:" << endl;

    for(int v=0; v<5; ++v){
        for(ll c=11; c<=20; ++c){
            string filename = file_path + "204101041_" + vowels[v] + "_" + to_string(c) + ".txt"; 
			vector<vector<ld>> file;
			file = readFile(filename);
			int fileSize = file.size();

			removeDCShift(file, fileSize);
			normalizeValues(file, fileSize);
			vector<vector<ld>> steadyFrames;
			steadyFrames = getSteadyFrames(file, fileSize);
			vector<vector<ld>> C;
			C = getCoefficients(steadyFrames);
			char predicted = predictVowel(C, codebook);
			if(predicted == vowels[v]) correct++;
			cout << predicted << "  ";
        }
		cout << endl;
	}

	cout << endl;
	cout << "Number of vowels predicted correctly: " << correct << "/50" << endl;
	cout << "Accuracy of the model: " << correct*1.0/50*100 << "%" << endl;
}

// Training the model on the training files and testing on the testing files
int main(){
	string file_path = "..\\Sound\\";
	vector<vector<ld>> codebook;
	cout << "Training the model..." << endl;
    codebook = train(file_path);
	cout << "Model trained." << endl;
	saveCodebook(codebook);
	printCodebook(codebook);
	cout << "\nTesting the model..." << endl;
	test(file_path, codebook);

	system("pause");
    return 0;
}