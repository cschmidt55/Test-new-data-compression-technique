%              f_Statistics_Set_Shaping_Theory_Huffman_c
%                        Contact info
%                Christian.Schmidt55u@gmail.com
%                  adrain.vdberg66@yahoo.com 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The program performs the following operations:
% 1) generates a random sequence with uniform distribution
% 2) calculate the frequencies of the symbols present in the sequence 
% 3) use this information to calculate the total information content sequence
% 4) apply the transform 
% 5) code the transformed sequence
% 6) compares the total information content of the generated sequence 
%    with the length of the encoded transfromated sequence
% 7) repeats all these steps a number of times defined by the parameter history
% 8) display the average values obtained
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Important
%
% If you change the length of the sequence and the number of the generated 
% symbols, you have to be careful that the huffman encoding approximates the 
% information content of the sequence by about one bit. if you take too 
% long sequences the huffman algorithm becomes very inefficient therefore, 
% it cannot detect the advantage obtained by the transformation.
% As a general rulen if you take ns symbols the length of the sequence 
% must be about 2*ns, in this case the Huffman encoding approximates 
% the information content of the sequence by about one bit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Definition
%
% The term information content Ic of a sequence is taken as defined in the SST.
% Example, given a sequence:
% 123445
%  
% We calculate the frequencies of the symbols present in the sequence.
%
%  symbol 1  frequencies 1/6
%  symbol 2  frequencies 1/6
%  symbol 3  frequencies 1/6
%  symbol 4  frequencies 2/6
%  symbol 5  frequencies 1/6
%
%  Ic=-log2(1/6)-log2(1/6)-log2(1/6)-log2(2/6)-log2(2/6)-log2(1/6)=13.5 bit
%
% 13.5 bit represents the minimum length of the coded sequence 
% in which the symbols have been replaced with uniquely decodable code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Results for different settings
%
%  The fSST2 and invfSST2 function only works when the number of symbols ns 
%  is greater or equal than 20 and less or equal than 500 and the length of 
%  the sequence N is greater than 40 and less than 1000. So, it is recommended 
%  to generate random sequences with a number of symbols and length greater than 
%  Ns=20 and N=40.
%
%  Ps=probability with which the transformed sequence can be encoded using
%  a uniquely decodable code (huffman coding) with lower bit number 
%  than the information content Ic of the initial sequence.
%  ns= number of symbols
%  N=lenght of the sequence
%
%  ns   N    Ps 
%  40   80   79%
%  50  100   83%
%  60  120   88%
% 500 1000   94%
%
%  As you can see, increasing the symbol number increases the probability Ps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
history=1000;
ns=40;
len=80;
cs=0;
totcodel=0;
totinfc=0;
tottinfc=0;
itnent=0;
infc=0;
tinfc=0;
itinfc=0;

for i=1:history

 % Generation of the sequence with a uniform distribution

 symbols=1:ns;
 prob(1,1:ns)=1/ns;
 seq=randsrc(1,len,[symbols; prob]);

 % Total information content sequence

 infc=0;

 for i2=1:len

  sy=seq(1,i2);
  fs=nnz(seq==sy)/len;
  infc=infc-log2(fs);

 end

 % Start trasformation

 mcodel=10000;

 nseq=fSST2(seq);

 % The new sequence is long nlen=len+1
 
 nlen=len+1;

 % total information content of the transformed sequence of length nlen

  tinfc=0;

  for i2=1:nlen

   sy=nseq(1,i2);
   fs=nnz(nseq==sy)/nlen;
   tinfc=tinfc-log2(fs);

  end 

 % Having transformed the sequence, we have to redefine the length of the vectors that are used in the encoding

 index=0;

 for i2=1:ns

  fs=nnz(nseq==i2)/nlen;

  if fs > 0

   index=index+1;    

  end

 end 

 c=zeros(index,1);
 vs=zeros(index,1);
 index=0;

 % We calculate the frequencies of the symbols in the transformed sequence

 for i2=1:ns

  fs=nnz(nseq==i2)/nlen;

  if fs > 0

   index=index+1;    
   c(index)=nnz(nseq==i2)/nlen;  
   vs(index)=i2;

  end

 end 

 % We code the transformed sequence

 counts=c;

 tdict=huffmandict(vs,counts);
 tcode=huffmanenco(nseq,tdict);

 bcode=de2bi(tcode);
 tcodel=numel(bcode);

 % If the length of the encoded message of the transformed sequence is less than the information content of the original sequence, we increase the counter cs by one

 if tcodel < infc

  cs=cs+1;

 end

 % We apply the inverse transform and we obtain the initial sequence

 iseq=invfSST2(nseq);

 % We check that the obtained sequence is equal to the initial sequence

 flag=isequal(seq,iseq);

 if flag == false

   fprintf('Error, sequence not equal to the initial sequence\n');

 end

 totcodel=totcodel+tcodel;
 totinfc=totinfc+infc;
 tottinfc=tottinfc+tinfc;

end

% We calculate the average of the information content of the generated sequences,the average of the information content of the transformed sequences and the average of the length of the encoded transformed sequence
  
 medinfc=totinfc/history;
 medcodel=totcodel/history;
 medtinfc=tottinfc/history;

% We calculate the percentage of sequences where the length of the encoded transformed sequence is less than the total information content of the generated sequence

 pcs=(cs/history)*100;

% We display the average values obtained

 fprintf('The average of the information content of the generated sequences\n');
 medinfc

 fprintf('The average of the length of the encoded transformed sequence\n');
 medcodel

 fprintf('The average of the information content of the transformed sequences\n');
 medtinfc

 fprintf('Number of sequences where the length of the encoded transformed sequence is less than the total information content of the generated sequence\n');
 cs

 fprintf('There is a percentage of %2.0f%% that length of the encoded transformed sequence is less than the total information content of the generated sequence\n',pcs);

