# Vowel-Recognition-By-preprocessing-of-raw-Speech-Data
 Vowel Recognition

Algorithm for vowel recognition:
Training Phase:
Step 1: Read the recording files:
	- The files are read framewise. Every frame is a collection of 320 samples

Step 2: Remove DC shift from the recordings: 
	- The DC shift values is precomputed on a seperate silence recording named "Silence.txt". 
	- This DC shift value is used for any new train/test recording

Step 3: Normalize the values:
	- Calculate the scaling factor by the formula: Scaling factor = 10000 * maximum sample value
	- Multiply every sample of the file with scaling factor

Step 4: Get steady framesof the recording:
	- Find the steady frames of the recording by analyzing the STE values for the frame
	- Steady frames are the frames which generally have high STE values
	- 5 steady frames are used in the program. 2 frames before and after the maximum STE frame are also considered as steady frames.

Step 5: Perform linear predictive analysis on the steady frames to get the linear predctive coefficients:
	- Apply hamming window on the frame
	- Calculate the auto-correlation terms for the frame
	- Apply Durbin's algorithm to calculate the linear predictive coefficients

Step 6: Calculate the Cepstral coefficients by linear transformation of the linear predicitve coefficients:
	- Multiply every linear predicitve coefficient of the frame by -1
	- Calculate the Cepstral coefficients by transformation of the linear predictive coefficients(formulaes are available in the notes)
	- Apply raised sine window on the cepstral coefficients

Step 7: Create the codebook using the previously calculated cepstral coefficients:
	- For every vowel 10 training files are available
	- Calculate the cepstral coefficients for every steady frame by the above described procedure
	- Calculte average of these cepstral coefficients framewise i.e. average of frame 1 of all the 10 training files of the vowel and so on
	- These framewise averages are the code vectors for the codebook
	- There are total 25 code vectors in the codebook. Each vowel correspoding to 5 code vectors
	- These code vectors are stored in the file "codebook.txt"

Testing phase:
Step 1: Calculate the cepstral coeffcients for every steady frame of every test file of the frame using the preocedure described above. 
	- There are total 50 test files available. Each vowel has 10 test files.

Step 2: Calculate the Tokhuras distance from code vectors in the codebook framewise:
	- Calculate the Tokhura distance framewise from code vectors of every vowel i.e. Tokhura distance of 1st steady frame from 1st code vector of vowel a, Tokhura distance of 1st steady frame from 1st code vector of vowel a and so on. Do this for every vowel code vector. 
	- Average these values of Tokhura distances i.e. Calculate average for Tokhura distance for 5 code vectors of vowel a, vowel e and so on.
	- The vowel which corresponds to the minimum average Tokhura's distance is the predicted vowwl.

Steps to execute the program:
	Step 1: Download and extract the file.
	Step 2: Open the project in Visual Studio 2010
	Step 3: To view the source code go to View -> Solution Explorer -> 204101041_VowelRecognition -> Source Files -> 204101041_VowelRecognition.cpp
	Step 4: Run the program by pressing F5 key. 
		The program will first train the model on the training files and create a code vector corresponding to the vowel recordings. This model will be then tested on the test recordings of every vowel. The program contains 10 training recordings for each vowel and 10 testing file for each vowel. The program creates a codebook of 25 code vectors using these 50 training files and this codebook is then used to predict the vowel in the test files. There are 10 testing files available for each vowel.  

Note:

1. While calculating the cepstral coefficients the linear predictive coefficients are not negate. To test the performance of the model by negating the ai uncomment the line number 155 of the code.
2. After calculating the cepstral coefficients raised sine window is not applied on the cepstral coefficients. To test the performance of the model by applying the raised sine window uncomment the line number 186 of the code.
3. To view the codebook open the file "codebook.txt" available in the project folder.
