The code of this paper:

T. Tang, C. Yang, S. Yan, L. Xu and D. Chen, "A Sparse Bayesian Learning-Based Approach With Indian Buffet Process Prior for Joint Wideband DOA and Frequency Band Estimation", IEEE Trans. Aerosp. Electron. Syst. .

%%%%%%

Description:

testAL1.m ： 		main function

SBLIBP.m :          function of our algorithm

%%%%%%

We provide the following simulated array received signals with two SNRs for algorithm verification.

1. receivedArrayData0dB.mat :  	Simulated array received data with 0dB SNR

2. receivedArrayData-5dB.mat : Simulated array received data with -5dB SNR

YUsed: array received data with size of  M* T* fft_num_use, 

where  

				M: the number of array elements.
                          
			    T:  the number of snapshots
                          
				fft_num_use:  the number of frequency points. The monitored frequency band is 400-2000 Hz                 

%%%%%%%

Signal parameters and array parameters have been shown in testAL1.m or the paper.  
