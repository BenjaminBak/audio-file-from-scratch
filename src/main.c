//
//  main.c
//  AudioFileFromScratch
//
//  Created by kado on 27.02.26.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define EPS 1e-12

const int sampleRate = 48000;
const int bitsPerSample = 16;
const int numChannels = 1;
const float durationInSeconds = 2.0;
const double frequency = 440.0;

typedef struct {
    char ChunkID[4];
    int ChunkSize;
    char Format[4];
    char SubChunk1ID[4];
    int SubChunk1Size;
    short AudioFormat;
    short NumChannels;
    int SampleRate;
    int ByteRate;
    short BlockAlign;
    short BitsPerSample;
    char SubChunk2ID[4];
    int SubChunk2Size;
} WavHeader;

void writeWavHeader(FILE* file, int nSamples) {
    WavHeader header = {
        /* RIFF/WAVE header */
        .ChunkID = "RIFF",
        .ChunkSize = 36 + nSamples * numChannels * bitsPerSample / 8,
        .Format = "WAVE",
        
        /* fmt subchunk */
        .SubChunk1ID = "fmt ",
        .SubChunk1Size = 16,
        .AudioFormat = 1,
        .NumChannels = numChannels,
        .SampleRate = sampleRate,
        .ByteRate = sampleRate * numChannels * bitsPerSample / 8,
        .BlockAlign = numChannels * bitsPerSample / 8,
        .BitsPerSample = bitsPerSample,
        
        /* data subchunk */
        .SubChunk2ID = "data",
        .SubChunk2Size = nSamples * numChannels * bitsPerSample / 8
    };
    fwrite(&header, sizeof(header), 1, file);
}

void writeWavData(FILE* file, double* samples, int numSamples) {
    short sample;
    for(int i = 0; i < numSamples; i++) {
        sample = (short)(samples[i] * (pow(2, bitsPerSample - 1) - 1));
        fwrite(&sample, sizeof(sample), 1, file);
    }
}

void writeSinWave(double* samples, int numSamples, float freq, float maxAmplitude) {
    double increment = 2.0 * M_PI * freq * 1.0/sampleRate;
    double phase = 0.0;
    for (int i = 0; i<numSamples; i++) {
        samples[i] = sin(phase)*maxAmplitude;
        phase += increment;
    }
}
void writeWhiteNoise(double* samples, int numSamples, float maxAmplitude) {
    for (int i = 0; i<numSamples; i++) {
        samples[i] = (((double)rand()/RAND_MAX)*2.0-1.0)*maxAmplitude;
    }
}

int applyEnv(double* samples, size_t numSamples, double* timeValPairs, size_t numTimeValPairs) {
    // timeValPairs[0] = pairOneTime; timeValPairs[0] = pairOneAmplitude;
    if (numTimeValPairs<1){
        return 0;
    }
    // check if pairs are in right order
    double prevTime = -INFINITY;
    for (size_t i = 0; i<numTimeValPairs; i++) {
        double t = timeValPairs[i*2];
        if ( t + EPS <prevTime) {
            fprintf(stderr,"Evenlope Pair in wrong order.\n");
            return 1;
        }
        prevTime = t;
    }
    size_t currentPair = 0;
    int state = 0; // 0 -> before first Value, 1 -> between values, 2 -> behind last value
    double timeA = timeValPairs[0];
    double valueA = timeValPairs[1];
    double timeB = timeValPairs[0];
    double valueB = timeValPairs[1];
    double valueDelta = 0.0;
    double timeDelta = 1.0;
    double interpolatedVal = 0.0;
    
    for (int i = 0; i<numSamples; i++) {
        double cTime = (double)i/(double)sampleRate;
        switch (state) {
            case 0:
                samples[i] *= valueA;
                if (timeA<=cTime) {
                    state = 1;
                    printf("0->1 at %f\n", cTime);
                } else break;
            case 1:
                if (cTime + EPS >= timeB) {
                    printf(" hmm ");
                    if (currentPair == numTimeValPairs - 1) {
                        state = 2;
                        interpolatedVal = valueB;
                        printf("1->2 at %f\n", cTime);
                    } else {
                        printf("switch from %lu to %lu\n",currentPair,currentPair+1);
                        
                        timeA = timeValPairs[currentPair*2];
                        valueA = timeValPairs[currentPair*2+1];
                        currentPair++;
                        timeB = timeValPairs[(currentPair)*2];
                        valueB = timeValPairs[(currentPair)*2+1];
                        valueDelta = valueB-valueA;
                        if (fabs(timeDelta) < EPS) {
                            timeDelta = 0.0;
                        }
                        timeDelta = timeB-timeA;
                    }
                } else {
                    double timeInSeg = cTime - timeA;
                    if (fabs(timeDelta) < EPS) {
                        interpolatedVal = valueB;
                    } else {
                        double t = timeInSeg / timeDelta;
                        interpolatedVal = valueA + valueDelta * t;
                    }
                }
            case 2:
                samples[i] *= interpolatedVal;
                break;
        }
        if (interpolatedVal>1.0) printf("!");
    }
    printf("px\n");
    return 0;
}



int main(int argc, const char * argv[]) {
    const int nSamples = sampleRate*durationInSeconds;
    double *samples = malloc(sizeof(double)*nSamples);
    writeWhiteNoise(samples, nSamples, 0.1);
    double env[6] = {0.0,0.0, 1.75,1.0, 2.0,0.0};
    if (applyEnv(samples, nSamples, env, 3)) return 1;
    
    FILE *fptr;
    fptr = fopen("monoWhiteNoiseWithEnv.wav", "wb");
    if (!fptr) {
        printf("Error opening file\n");
        return 1;
    }
    writeWavHeader(fptr, nSamples);
    writeWavData(fptr, samples, nSamples);
    fclose(fptr);
    free(samples);
    return EXIT_SUCCESS;
}
