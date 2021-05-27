#include <iostream>
#include <stdio.h>
#include <fstream>
#include <list>
#include <iomanip>
#include <sstream>

//Matrix stuff
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

using namespace std; 

void checkString();
void parseInput(string);
float stringToFloat(string);
void initialize(string);

//Structs that should probably be in a header...
struct sphere{
	string name;
	glm::vec3 pos;
	glm::vec3 scl;
	glm::vec3 color;
	float ka, kd, ks, kr, n;
};

struct light{
	string name;
	glm::vec3 pos;
	glm::vec3 color;
};

struct ray{
	float t;
	glm::vec3 s;
	glm::vec3 c;
	
};

//The sphere and light lists
list <sphere> spheres;
list <light> lights;

//All the values that will be parsed from the text file
float NEAR;
float LEFT;
float RIGHT;
float BOTTOM;
float TOP;
int RES[2];
float BACK[3];
glm::vec3 AMBIENT;
string OUTPUT;

//The number of reflects
int reflects=0;

glm::vec3 eye = glm::vec3(0,0,0);
glm::vec4 u = glm::vec4(1,0,0,0);
glm::vec4 v = glm::vec4(0,1,0,0);
glm::vec4 n = glm::vec4(0,0,1,0);

//the current line to be parsed
string parseArray[16];

//Used as a setter for the sphere struct
//Returns the sphere to be added to a list
sphere setSphere(){
	struct sphere sph;
	
	sph.name = parseArray[1];
	sph.pos = glm::vec3(stof(parseArray[2]), stof(parseArray[3]), stof(parseArray[4]));
	sph.scl = glm::vec3(stof(parseArray[5]), stof(parseArray[6]), stof(parseArray[7]));
	sph.color = glm::vec3(stof(parseArray[8]), stof(parseArray[9]), stof(parseArray[10]));
	sph.ka = stof(parseArray[11]); sph.kd = stof(parseArray[12]); sph.ks = stof(parseArray[13]); sph.kr = stof(parseArray[14]); sph.n = stof(parseArray[15]);
	
	return sph;
}

//Used as a setter for the light struct
//Returns the light to be added to a list
light setLight(){
	struct light lig;
	
	lig.name = parseArray[1];
	lig.pos = glm::vec3(stof(parseArray[2]), stof(parseArray[3]), stof(parseArray[4]));
	lig.color = glm::vec3(stof(parseArray[5]), stof(parseArray[6]), stof(parseArray[7]));
	return lig;
}

//Creates a ray where eye is the start, c is the vector,
//and t is the distance and returns it w
ray raycast(glm::vec3 eye, glm::vec3 c, float t){
	//return eye + (t *(eye + p));
	struct ray r;
	r.t = t;
	r.c = c;
	r.s = eye;
	
	return r;
}

//A gross parser function that checks the current input array
void checkString(){
	
	if(parseArray[0].compare("NEAR") == 0){
		NEAR = stof(parseArray[1]);

	}
	else if(parseArray[0].compare("LEFT") == 0){
		LEFT = stof(parseArray[1]);
	}
	else if(parseArray[0].compare("RIGHT") == 0){
		RIGHT = stof(parseArray[1]);
	}
	else if(parseArray[0].compare("BOTTOM") == 0){
		BOTTOM = stof(parseArray[1]);
	}
	else if(parseArray[0].compare("TOP") == 0){
		TOP = stof(parseArray[1]);
	}
	else if(parseArray[0].compare("RES") == 0){
		RES[0] = stoi(parseArray[1]);
		RES[1] = stoi(parseArray[2]);
	}
	else if(parseArray[0].compare("SPHERE") == 0){
		spheres.push_back(setSphere());
	}
	else if(parseArray[0].compare("LIGHT") == 0){
		lights.push_back(setLight());
	}
	else if(parseArray[0].compare("BACK") == 0){
		BACK[0] = stof(parseArray[1]);
		BACK[1] = stof(parseArray[2]);
		BACK[2] = stof(parseArray[3]);
	}
	else if(parseArray[0].compare("AMBIENT") == 0){
		AMBIENT = glm::vec3(stof(parseArray[1]),stof(parseArray[2]),stof(parseArray[3]));
	}
	else if(parseArray[0].compare("OUTPUT") == 0){
		OUTPUT = parseArray[1];
	}
}

//Ray intersection function that returns a vec2 of both hits
glm::vec2 sphereIntersection(ray r){
	
	float A = glm::dot(r.c, r.c);
	float B = glm::dot(r.s, r.c);
	float C = glm::dot(r.s, r.s)-1;
	
	float temp = ((B*B) - A*C);
	if(temp<0){
		return glm::vec2(-1, -1);
	}
	temp = (sqrt(temp)/A);
	
	double inter1 = (-(B/A)) + temp;
	double inter2 = (-(B/A)) - temp;
	
	if(inter1<inter2){
		return glm::vec2(inter1, inter2);
	}
	return glm::vec2(inter2, inter1);
}

bool inShadow(glm::vec3 rayHit, light lig, sphere currSphere, glm::mat4 tempM, glm::vec3 cannHit){
	
	glm::mat4 M;
	for(const auto& s: spheres){
		if(currSphere.name.compare(s.name) == 0){
			//cout << "return\n";
			continue;
		}
		
		//The current sphere intersections matrix
		//I could probably save computation by having this stored in the struct
		M = glm::mat4(1.0);
		M = glm::translate(M, s.pos);
		M = glm::scale(M, s.scl);
		M = glm::inverse(M);
		
		//Creates the shadow vector
		glm::vec3 shadowVec = (lig.pos-rayHit);
		shadowVec = glm::normalize(shadowVec);
	
		//creating the shadow ray
		ray shadowRay = raycast(rayHit, shadowVec, 80);
		
		//The s and c position of hit
		glm::vec4 hS = glm::vec4(rayHit,1.0);
		glm::vec4 hC = glm::vec4(shadowVec,0.0);
		
		//transforming the hit
		hS = M*(hS);
		hC = M*shadowRay.t*(hC);
		
		shadowRay.s = hS;
		shadowRay.c = hC;
		
		//checks the vector
		float temp = sphereIntersection(shadowRay)[0];
	
		//Makes sure the sphere is in shadow with an epsilon for float error
		float epsilon = 0.00001;
		if(temp>=0-epsilon){
			return true;
		}
	}
	return false;
}

//The raytracing function that can be called recursively
glm::vec3 raytrace(ray r){
	glm::vec3 color = glm::vec3(0,0,0);
	
	//Transform matrix
	glm::mat4 M;
	glm::mat4 iM;
	
	//Normal
	glm::vec4 N = glm::vec4(0,0,0,1.0);
	glm::vec3 L = glm::vec3(0,0,0);
	glm::vec3 R = glm::vec3(0,0,0);
	glm::vec3 V = glm::vec3(0,0,0);	
	
	//A hard coded big number to use as the back
	float closestIntersect=10000;

	//If its the first reflect, print the background color
	if(reflects==0){
		color = glm::vec3(BACK[0],BACK[1],BACK[2]);
	}
	//How many reflections will happen 
	if(reflects >=3){
		return glm::vec3(0,0,0);
	}
	
	for(const auto& s: spheres){
		
		//Used to store the past r.c to revert after transform
		glm::vec3 tempC = r.c;
		glm::vec3 tempS = r.s;
				
		//transformation matrix
		M = glm::mat4(1.0);
		M = glm::translate(M, s.pos);
		M = glm::scale(M, s.scl);
		iM = glm::inverse(M);
		
		//homogeneous S and C
		glm::vec4 hS = glm::vec4(r.s,1.0);
		glm::vec4 hC = glm::vec4(r.c,0.0);
		
		//Transforming by the inverse matrix
		hS = (iM * hS);		
		hC = r.t*(iM * (hC));
		
		//converting from homogeneous and normalizing
		r.s = glm::vec3(hS);
		r.c = glm::normalize(glm::vec3(hC));
		
		//Check the intersection
		glm::vec2 vecTh = sphereIntersection(r);
		float th = vecTh[0];
		
		if(th>=(0)){
			
			hS = glm::vec4(r.s,1.0);
			hC = glm::vec4(r.c,0.0);
			
			//The rayhit on the transformed sphere and cannonical sphere
			glm::vec4 rayHit = M*hS+(M*th*hC);
			glm::vec3 rayHit3 = glm::vec3(rayHit);
			glm::vec4 cannHit = hS+(th*hC);
			
			//Check the intersection
			if(glm::length(rayHit3) < closestIntersect){
				
				//Used for the checking the closest intersection
				closestIntersect = glm::length(rayHit3);
			
				//Subtracting the origin doesnt do anything, but just incase I change stuff later
				N = (cannHit - glm::vec4(0, 0, 0, 0));
				N = glm::transpose(iM)*N;
				glm::vec3 uNorm = glm::vec3(N);
				uNorm = glm::normalize(uNorm);
			
				//Ambient
				color = s.ka * AMBIENT * s.color;
				
				//For all lights
				for(const auto& l: lights){	
				
					//If in shadow
					if(!inShadow(rayHit, l, s, iM, cannHit)){
						
						//Vectors used for diffuse and specular
						L = l.pos-rayHit3;
						L = glm::normalize(L);
						R = glm::reflect(-L, uNorm);
						R = glm::normalize(R);
						V = eye-rayHit3;
						V = glm::normalize(V);
						
						//Diffuse
						color += glm::clamp(s.kd * l.color * (glm::dot(uNorm,L)* s.color),0.0f, 1.0f);
					
						//Specular
						if(glm::dot(L,uNorm) >= 0.0f){
							color += glm::clamp((s.ks * l.color * pow(max(glm::dot(R,V),0.0f),s.n)),0.0f, 1.0f);
						}
					}
				}
				
				//Reflections
				if(s.kr > 0.1){
					glm::vec3 vrf = -2.0f*(glm::dot(uNorm,tempC)*uNorm)+tempC;
					ray reflRay;
					reflects++;
					reflRay = raycast(rayHit3, vrf, 80);
					color += glm::clamp(s.kr * raytrace(reflRay),0.0f,1.0f);
				}
			}
		}
		//undo transformations
		r.c = tempC;
		r.s = tempS;
	}
	color = glm::clamp(color, 0.0f, 1.0f);
	return color;
}

//Basic initialization and string parsing to be passed to checkString
void initialize(string inputFile){
	ifstream file(inputFile);
	string str;

	int i=0;
	
	//Splits the file into a string of the current line
	while(std::getline(file, str)){
		istringstream iss(str);
		//If a line is empty continue
		if(str.compare("") == 0){
			continue;
		}
		
		//parse the line into words and place them in an array
		string word;
		while(iss >> word){
			parseArray[i] = word;
			i++;
		}
		i=0;
		
		//Convert the array into the global variables
		checkString();
	}
}

//Do the math for the pixel coordinates
float pixelCoord(int nSize, int rc, float p, float n){
	return n+(p*(2*rc)/nSize);
}

// Output in P6 format, a binary file containing:
// P6
// ncolumns nrows
// Max colour value
// colours in binary format thus unreadable
void save_imageP6(int Width, int Height, string OUTPUT, unsigned char* pixels) {
	char* fname;
	fname = &OUTPUT[0];
	
	FILE *fp;
	const int maxVal=255; 
  
	printf("Saving image %s: %d x %d\n", fname,Width,Height);
	fp = fopen(fname,"wb");
	if (!fp) {
        printf("Unable to open file '%s'\n",fname);
        return;
	}
	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", Width, Height);
	fprintf(fp, "%d\n", maxVal);

	for(int j = 0; j < Height; j++) {
		fwrite(&pixels[j*Width*3], 3,Width,fp);
	}
	fclose(fp);
}

int main(int argc, char *argv[]){
	//Checks to see if the program is running and initilializes it
	std::cout << "Running\n";
	if (argc < 2) {
		cout << "No input file exiting\n";
		exit(1);
	}
	
	initialize(argv[1]);

	
	unsigned char *pixels;
	pixels = new unsigned char [3*RES[0]*RES[1]];
	//Image plane in camera coords
	glm::vec3 imagePl = glm::vec3(0,0,-NEAR);
	
	//Transform matrix
	glm::mat4 M(1.0);
	//t value
	float t = 80;
	//The ray struct to be used
	struct ray r;
	
	int k=0;
	//For every pixel on the screen
	for(int i=0; i<RES[1]; i++){
		//vr
		imagePl[1] = pixelCoord(RES[1], i, LEFT, RIGHT);
		for(int j=0; j<RES[0]; j++){
			//uc
			imagePl[0] = pixelCoord(RES[0], j, TOP, BOTTOM);

			r = raycast(eye, imagePl, t);

			//Used to color the pixels in the image	
			glm::vec3 color = raytrace(r);
			reflects = 0;
			
			pixels[k] = color[0]*255;
			pixels[k+1] = color[1]*255;
			pixels[k+2] = color[2]*255;
				
			k+=3;
		}
		cout << "here ";
	}
	
	//Saves the image to the file OUTPUT
	save_imageP6(RES[0], RES[1], OUTPUT, pixels);
	return 0;
}
