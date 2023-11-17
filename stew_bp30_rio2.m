% expandable no swap!!!

%%%%%%%%
%%%%%%%%
%%%%%%%%
function stew_bp30();
global blkcountINSRT wmCount; % INSERT
global BstreamBlockCount BSTREAM  RIFCAEP_CNT CarrayCNT; % RECOVER
% INSERT
blkcountINSRT=1;
wmCount=1;
% RECOVER
BstreamBlockCount=1;
BSTREAM=[];
RIFCAEP_CNT=1;
CarrayCNT=1;

% MAKES HUGE DIFFERENCE IF USE ReadImage Lenna or imread lena.jpg
lena=ReadImage('Lenna'); 
%lena=MyReadImage('Lenna'); 
%lena=double(imread('lena_gray_256.jpg'));
%lena=lena+1;
%lena=ReadImage('Daubechies');
%lena=lena(101:164,101:164);
%lena=lena(1:10,1:4);
%lena=[lena; 255 255 0 0;20 20 106 100;100 106 21 20;20 21 21 21;105 100 1 0;0 2 0 2;255 253 255 253]; 

% NOTE: Lenna256: th1=17 and th2=27, bpp9=0.09+0.34=0.43, psnr=38.82: same results as Tian for lena512
% I am doing horizontal pair 1st (weakest) then doing the vertical pairs (stronger)
% Doing weakest then strongest give best results (all the time???)
 th1=21; th2=35;
%th1=30; th2=45;

% create a payload
PL=[0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
PL=[PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL];
PL=[PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL];
PL=[PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL];
PL=[PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL PL];
PL=randn(size(PL)) > 0;

[lenaX,bi1,pi1,bpp1]=INSERT21(lena,th2,PL);
blkcountINSRT=1; wmCount=1;  % reset global variables

	%[lenaXX,bi2,pi2,bpp2]=INSERT21(lenaX,th2,PL);
	%blkcountINSRT=1; wmCount=1;  % reset global variables
	%[LX,br2,pr2]=RECOVER21(lenaXX);
	%BstreamBlockCount=1; BSTREAM=[]; RIFCAEP_CNT=1; CarrayCNT=1; % reset globals

[LR,br1,pr1]=RECOVER21(lenaX);
BstreamBlockCount=1; BSTREAM=[]; RIFCAEP_CNT=1; CarrayCNT=1; % reset globals

% stats
theerror=lena - LR;
fprintf('\nlena and LR (Lena Recovered): max error=%.2f\n',max(max(theerror)));
fprintf('lena and LR (Lena Recovered): min error=%.2f\n',min(min(theerror)));

fprintf('\nNC between inserted B Streams [bi1] and recovered B streams [br1]\n'); 
print_NC(bi1,br1);
fprintf('length BSTREAM [bi1] = %d\n',length(bi1));
fprintf('NC between insertedPayload [pi1] and recovered payload [pr1]\n'); 
print_NC(pi1,pr1);
fprintf('length payload [pi1] = %d\n',length(pi1));

figure, GrayImage([lena lenaX LR]);
title('Original    WatermarkedOnce   Recovered');

%fprintf('\nNC between inserted B Streams [bi2] and recovered B streams [br2]\n'); 
%print_NC(bi2,br2);
%fprintf('length BSTREAM [bi2] = %d\n',length(bi2));
%fprintf('NC between insertedPayload [pi2] and recovered payload [pr2]\n'); 
%print_NC(pi2,pr2);
%fprintf('length payload [pi2] = %d\n',length(pi2));

%figure, GrayImage([lena lenaX LR]);
%title('Original    WatermarkedTwice    Recovered');
fprintf('\npsnr=%.2f\n',psnrMetric(lena,lenaX,2));
%fprintf('total bpp = bpp1 + bpp2 = %.2f + %.2f = %.2f\n',bpp1,bpp2,bpp1+bpp2);

fprintf('\nth1=%.2f th2=%.2f\n',th1,th2);

fprintf('PayLoad [PL] 1 to 10 (checking randomness!)\n');
for i=1:10,
	fprintf('%d',PL(i));
end
fprintf('\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lenaDashEC,insertedStream,insertedPayload,bpp9]=INSERT12(lena,th,PL);
global blkcountINSRT wmCount;
% notes:
% h=-2, 253 255 and 0 2, changeable, rest expandable! test with wm=0 and 1
% 255 255 pair, not expandable or changeable.

% classify
expandablePairMatrix=blkproc(lena,[1 2],@classifyExpandableFun,th);
changeablePairMatrix=blkproc(lena,[1 2],@classifyChangeableFun);
expandablePairArray=expandablePairMatrix';
expandablePairArray=ShapeAsRow(expandablePairArray(:));
changeablePairMatrix=changeablePairMatrix & ~expandablePairMatrix;
changeablePairArray=changeablePairMatrix';
changeablePairArray=ShapeAsRow(changeablePairArray(:));

totalNoOfCEpairs = sum(changeablePairArray) + sum(expandablePairArray);

hValueMatrix=blkproc(lena,[1 2],@getThe_h_DifferenceValuesFun);
Carray=getCarray(hValueMatrix,changeablePairMatrix); % store the original LSBs of changed pairs

imwrite(uint8(expandablePairArray*255),'exLM.pbm');
!pbmtojbg exLM.pbm exLM.jbg
fp=fopen('exLM.jbg','rb');
arrSymbols=fread(fp,'uchar');
fprintf('\nlength of location map=%d BYTES, %d*8=%d BITS %d*9=%d BITS\n', ...
	length(arrSymbols),length(arrSymbols),length(arrSymbols)*8,length(arrSymbols),length(arrSymbols)*9);
fclose(fp);
len9=[];
for i=1:length(arrSymbols),
	theValArr=dec2bin(arrSymbols(i),9);
	nineBin=double(uint8(theValArr)-48);
	len9=[len9 nineBin];
end
EOS=[1 1 1 1 1 1 1 1 1];
LarrComp=[len9 EOS]';  % locmap compressed as binary stream



% arrSymbols is compressed location map (theoretical just now!)
LC8 = length(arrSymbols)*8 + 8 + length(Carray); % amount of pairs needed for compressed LOCATION MAP and Carray (8bits)
LC9 = length(arrSymbols)*9 + 9 + length(Carray); % amount of pairs needed for compressed LOCATION MAP and Carray (9bits)
plBITS8 = totalNoOfCEpairs - LC8; % amount of pairs for embedding payload (no EOS, 8 bits)
plBITS9 = totalNoOfCEpairs - LC9; % amount of pairs for embedding payload (with EOS, 9 bits)
bpp8=plBITS8/prod(size(lena));
bpp9=plBITS9/prod(size(lena));
fprintf('\nbpp8=%.2f bpp9=%.2f\n',bpp8,bpp9);
fprintf('\nnoOfExpandablePairs=%d, numberOfChangeablePairs=%d\n',sum(expandablePairArray),sum(changeablePairArray));
fprintf('lengthPayload [8bits, no EOS]=%d, lengthPayload [9bits, with EOS]=%d\n',plBITS8,plBITS9);
fprintf('totalNoOfPairs=%d, TotalNoOfExpandableAndChangeablePairs=%d CarrayLength=%d\n',prod(size(lena))/2,totalNoOfCEpairs,length(Carray));
fprintf('length compressed location map = %d bits\n',length(LarrComp));
locmapsize=prod(size(lena))/2;

% concatenate the watermark ... LCP
wm=[LarrComp;Carray;PL(:)];

% insert changeable and expandable pairs
lenaDashEC=blkproc(lena,[1 2],@insertChangeableAndExpandable,changeablePairArray,expandablePairArray,wm);

fprintf('\nmax lena=%.2f  min lena=%.2f\n',max(lena(:)),min(lena(:)));
fprintf('max lenaDashEC=%.2f  min lenaDashEC=%.2f\n\n',max(lenaDashEC(:)),min(lenaDashEC(:)));

insertedStream = wm(1:totalNoOfCEpairs); % stream inserted into LSBs of changeable and expandable
%insertedPayload=insertedStream(length(Carray)+1:length(insertedStream));
lengthPL=totalNoOfCEpairs-length(LarrComp)-length(Carray);
insertedPayload=wm(length(LarrComp)+length(Carray)+1:totalNoOfCEpairs);

figure, GrayImage(expandablePairMatrix*255), title('Location map');



%%%%%%%%%%%%%%%%%%%% RECOVER %%%%%%%%%%%%%%%%%%%
% classify
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
function [LR,theBstreamRec,recoveredPayload]=RECOVER12(lenaDashEC);
global BstreamBlockCount BSTREAM  RIFCAEP_CNT CarrayCNT;
changeablePairMatrix2=blkproc(lenaDashEC,[1 2],@RECOVERclassifyChangeableFun);
changeablePairArray2=changeablePairMatrix2';
changeablePairArray2=ShapeAsRow(changeablePairArray2(:));

% FIND ALL CHANGEABLE PAIRS, store as B stream (BSTREAM)
dummy=blkproc(lenaDashEC,[1 2],@B_STREAM_recoverEC_LSBs,changeablePairArray2);
theBstreamRec=BSTREAM;

fprintf('length theBstreamRec=%d\n',length(theBstreamRec));

% expandable pairs are changeable!!
% would need to get the compressed locmap from B stream, decompress it, then 
% i have a map of where all the expanded pairs are. Here, it is just copied,
% assumed i have decompressed it
%LOCMAP=expandablePairMatrix;
%LOCMAP=LOCMAP';
%LOCMAP=ShapeAsRow(LOCMAP(:));
% see as.m and as2.m
asr=[];
thetest=0;
i=1;
%for i=1:9:length(len9),
while thetest == 0,
	get9=BSTREAM(i:i+8);
	kk=strread(num2str(get9),'%s')';
	kk=char(kk)';
	decval=double(bin2dec(kk));
	%fprintf('i=%d  decval=%d\n',i,decval);
	if decval == 511
		fprintf('EOS found\n');
		BSTREAM=BSTREAM(i+9:length(BSTREAM)); % chop BSTREAM (remove L, CP left)
		thetest=1;
		break;
	end
	asr=[asr decval];
	i=i+9;
end
fp2=fopen('lmrec.jbg','wb');
fwrite(fp2,asr,'uchar');
fclose(fp2);
!jbgtopbm lmrec.jbg lmrec.pbm
theLM=imread('lmrec.pbm');
LOCMAP=double(theLM);
fprintf('length BSTREAM=%d\n',length(BSTREAM));
fprintf('length LOCMAP=%d\n',length(LOCMAP));

% from all the changeable pairs, pick out the exapandable pairs
% changeablePairArray2 contains locations of all changeable and
% expandable pairs (expandable pairs are changeable!)
recoverExpandedArray = changeablePairArray2.*LOCMAP;
% pick out the pairs that were changed
recoverChangedArray=changeablePairArray2 .* ~LOCMAP;

% recover the original pair values that were changed and the original pair values that were expanded.
LR=blkproc(lenaDashEC,[1 2],@recoverImageFromChangeableAndExpandablePairs,recoverChangedArray,recoverExpandedArray,BSTREAM);
% stats
% CarrayCNT points to one place beyond length of Carray (see recoverImageFromChangeableAndExpandablePairs function).
%recoveredPayload=theBstreamRec(CarrayCNT:length(theBstreamRec)); 
recoveredPayload=BSTREAM(CarrayCNT:length(BSTREAM)); 



%%%%%%
%%%%%%
function blkout=recoverImageFromChangeableAndExpandablePairs(blkin,changeablePairArray,recoverExpandedArray,theBstr);
global RIFCAEP_CNT CarrayCNT;
blkout=zeros(size(blkin));
if changeablePairArray(RIFCAEP_CNT)==1,
	x=blkin(1);
	y=blkin(2);
	%
	% NO SWAP!
	% 
	xx=x; yy=y;
	l = floor((xx+yy)/2);
	h = xx - yy;
	% have found a pair that was changed
	if h==0 | h==1,
		hnew=1;
	elseif h==-2 | h==-1
		hnew=-2;
	else
		thebit=theBstr(CarrayCNT);
		CarrayCNT=CarrayCNT+1;
		hnew=(2*floor(h/2))+thebit;
	end
	xdash=l+floor((hnew+1)/2);
	ydash=l-floor(hnew/2);
	%
	% NO SWAP
	%
	%blkout=[xdash ydash];
	blkout(1)=xdash;
	blkout(2)=ydash;
elseif recoverExpandedArray(RIFCAEP_CNT)==1,
	% found an expanded pair
	x=blkin(1);
	y=blkin(2);
	%[maxval,maxloc]=max([x y]);
	%swap=0;
	%if maxloc == 1,
	%	swap=0;
	%	xx = x;
	%	yy = y;
	%elseif maxloc == 2,
	%	swap = 1;
	%	xx = y;
	%	yy = x;
	%end % if swap is 1, then orginal h value was NEGATIVE
	xx=x; yy=y;
	l = floor((xx+yy)/2);
	h = xx - yy; % NEW
	% now decode to get original image back
	hnew=floor(h/2); % rule for getting hnew in EXPANDED pair
	xdash=l+floor((hnew+1)/2);
	ydash=l-floor(hnew/2);
	%if swap == 1,	% original h value of this pair was NEGATIVE
	%	%blkout=[ydash xdash];
	%	blkout(1)=ydash;
	%	blkout(2)=xdash;
	%else
	%	%blkout=[xdash ydash];
	%	blkout(1)=xdash;
	%	blkout(2)=ydash;
	%end
	blkout(1)=xdash;
	blkout(2)=ydash;
	%fprintf('x=%d y=%d l=%d h=%d hnew=%d xdash=%d ydash=%d swap=%d\n',x,y,l,h,hnew,xdash,ydash,swap);
else
	% do not change nor expand the pair
	blkout=blkin;
end
RIFCAEP_CNT=RIFCAEP_CNT+1;

%%%%%%
%%%%%%
function Carray=getCarray(hValueMatrix,changeablePairMatrix);
Cmat=hValueMatrix .* changeablePairMatrix;
Cmat=Cmat';
Cmat=ShapeAsRow(Cmat(:));
Carray=[];
% remove 1 and -2
for i=1:length(Cmat),
	%fprintf('cmat i = %d\n',Cmat(i));
	if (Cmat(i) ~= 1) & (Cmat(i)~=-2),
		Carray=[Carray; Cmat(i)];
	end
end
Carray=nonzeros(Carray);
Carray=mod(abs(Carray),2);

%%%%%
%%%%%
function blkout = RECOVERclassifyChangeableFun(blkin);
x=blkin(1);
y=blkin(2);
%[maxval,maxloc]=max([x y]);
%swap=0;
%if maxloc == 1,
%	swap=0;
%	xx = x;
%	yy = y;
%elseif maxloc == 2,
%	swap = 1;
%	xx = y;
%	yy = x;
%end
%
% NO SWAP!!!!!!!!!!!!!!!!!!!!!!!!
%
% floor(-5/2)=floor(-2.5), round down to -3!!
%
xx=x;yy=y;
l = floor((xx+yy)/2);
h = xx - yy;
% classify l value, changeable? 
LHS=2*(255-l);
RHS=(2*l)+1;
MINVAL_LEFT_RIGHT=min(LHS,RHS);
% changeable check
test0=abs((2*floor(h/2))+0) <= MINVAL_LEFT_RIGHT;
test1=abs((2*floor(h/2))+1) <= MINVAL_LEFT_RIGHT;
if (test0 & test1),
	outval = 1;
else
	% pair not changeable
	outval = 0;
end
blkout=outval;
if outval ~= 1
	%fprintf('RCC xx=%d yy=%d l=%d h=%d LHS=%d RHS=%d test0=%d test1=%d MINVAL=%d outval=%d\n',xx,yy,l,h,LHS,RHS,test0,test1,MINVAL_LEFT_RIGHT,outval);
end

%------------------
%------------------
% collect all LSBs from difference numbers h
function blkout=B_STREAM_recoverEC_LSBs(blkin,changedPairsArray);
global BstreamBlockCount BSTREAM;
blkout=[];
if changedPairsArray(BstreamBlockCount)==1,
	blkout=blkin;
	% found a LSB therefore store copy of it in BSTREAM
	x=blkin(1);
	y=blkin(2);
	%[maxval,maxloc]=max([x y]);
	%swap=0;
	%if maxloc == 1,
	%	swap=0;
	%	xx = x;
	%	yy = y;
	%elseif maxloc == 2,
	%	swap = 1;
	%	xx = y;
	%	yy = x;
	%end % if swap is 1, then orginal h value was NEGATIVE
	xx=x; yy=y;
	l = floor((xx+yy)/2);
	h = abs(xx - yy); % NEW, just find the magnitude of diff, odd: wm=1, even: wm=0
	binval=dec2bin(h);
	binval=double(uint8(binval(length(binval)))-48); % need double so no warnings
	BSTREAM=[BSTREAM binval]; % belongs set {0,1}
else	
	% wasn't expanded
	%fprintf('BSTREAM Recovery: not expanded!!!!!  BstreamBlockCount=%d\n\n\n\n',BstreamBlockCount);
	blkout=blkin;
end
BstreamBlockCount=BstreamBlockCount+1;

%------------------
%------------------
function blkout=insertChangeableAndExpandable (blkin,changeablePairArray,expandablePairArray,wm);
global blkcountINSRT wmCount;
blkout=zeros(size(blkin));
%blkout=[];
if changeablePairArray(blkcountINSRT)==1,
	x=blkin(1);
	y=blkin(2);
	%
	% NO SWAP!
	%
	xx=x; yy=y;
	l = floor((xx+yy)/2);
	h = xx - yy;
	hnew=(floor(h/2)*2) + wm(wmCount);
	xdash=l+floor((hnew+1)/2);
	ydash=l-floor(hnew/2);
	%
	% NO SWAP!
	%
	%blkout=[xdash ydash];
	blkout(1)=xdash;
	blkout(2)=ydash;
	%if xdash <= 0 | ydash <= 0,
	%	fprintf('CHANGEABLE InsertChangeable: x=%d y=%d xx=%d yy=%d l=%d h=%d hnew=%d wmval=%d xdash=%d ydash=%d\n',...
	%		x,y,xx,yy,l,h,hnew,wm(wmCount),xdash,ydash);
	%end
	wmCount=wmCount+1;
elseif expandablePairArray(blkcountINSRT)==1,
	x=blkin(1);
	y=blkin(2);
	%[maxval,maxloc]=max([x y]);
	%swap=0;
	%if maxloc == 1,
	%	swap=0;
	%	xx = x;
	%	yy = y;
	%elseif maxloc == 2,
	%	swap = 1;
	%	xx = y;
	%	yy = x;
	%end
	xx=x;
	yy=y;
	l = floor((xx+yy)/2);
	h = xx - yy; 
	hnew= (2*h) + wm(wmCount);
	xdash=l+floor((hnew+1)/2);
	ydash=l-floor(hnew/2);
	%if swap == 1,
	%	%blkout=[ydash xdash];
	%	blkout(1)=ydash;
	%	blkout(2)=xdash;
	%else
	%	%blkout=[xdash ydash];
	%	blkout(1)=xdash;
	%	blkout(2)=ydash;
	%end
	blkout(1)=xdash;
	blkout(2)=ydash;
	%if xdash <= 0 | ydash <= 0,
	%	fprintf('EXPANDABLE InsertChangeable: x=%d y=%d xx=%d yy=%d l=%d h=%d hnew=%d wmval=%d xdash=%d ydash=%d\n', ...
	% 		x,y,xx,yy,l,h,hnew,wm(wmCount),xdash,ydash);
	%end
	wmCount=wmCount+1;
else
	blkout=blkin;
end

blkcountINSRT=blkcountINSRT+1;

%------------------
%------------------
function blkout=getThe_h_DifferenceValuesFun(blkin);
x=blkin(1);
y=blkin(2);
[maxval,maxloc]=max([x y]);
%swap=0;
%if maxloc == 1,
%	swap=0;
%	xx = x;
%	yy = y;
%elseif maxloc == 2,
%	swap = 1;
%	xx = y;
%	yy = x;
%end
xx=x; yy=y;
l = floor((xx+yy)/2);
h = xx - yy;
blkout=h;

%%%%%
%%%%%
function blkout = classifyChangeableFun(blkin);
x=blkin(1);
y=blkin(2);
%[maxval,maxloc]=max([x y]);
%swap=0;
%if maxloc == 1,
%	swap=0;
%	xx = x;
%	yy = y;
%elseif maxloc == 2,
%	swap = 1;
%	xx = y;
%	yy = x;
%end
%
% NO SWAP!!!!!!!!!!!!!!!!!!!!!!!!
%
% floor(-5/2)=floor(-2.5), round down to -3!!
%
xx=x;yy=y;
l = floor((xx+yy)/2);
h = xx - yy;
% classify l value, changeable? 
LHS=2*(255-l);
RHS=(2*l)+1;
MINVAL_LEFT_RIGHT=min(LHS,RHS);
% changeable check
test0=abs((2*floor(h/2))+0) <= MINVAL_LEFT_RIGHT;
test1=abs((2*floor(h/2))+1) <= MINVAL_LEFT_RIGHT;
if (test0 & test1),
	if h==0 | h==-1,
		outval=0;
	else
		outval = 1;
	end
else
	% pair not changeable
	outval = 0;
end
blkout=outval;

%------------------
%------------------
function blkout = classifyExpandableFun(blkin,th);
x=blkin(1);
y=blkin(2);
%[maxval,maxloc]=max([x y]);
%swap=0;
%if maxloc == 1,
%	swap=0;
%	xx = x;
%	yy = y;
%elseif maxloc == 2,
%	swap = 1;
%	xx = y;
%	yy = x;
%end
xx=x; yy=y;
l = floor((xx+yy)/2);
h = xx - yy;
% classify l value, expandable?
LHS=2*(255-l);
RHS=(2*l)+1;
MINVAL_LEFT_RIGHT=min(LHS,RHS);
% expandable check
test0=abs((2*h)+0) <= MINVAL_LEFT_RIGHT;
test1=abs((2*h)+1) <= MINVAL_LEFT_RIGHT;
if test0 & test1,
	if abs(h) <= th,
		% pair are expandable
		outval=1;
	else
		% pair difference greater than th, do not expand
		outval = 0;
	end
else
	% pair not expandable
	outval = 0;
end
%fprintf('xx=%d yy=%d l=%d h=%d LHS=%d RHS=%d test0=%d test1=%d MINVAL=%d outval=%d\n',xx,yy,l,h,LHS,RHS,test0,test1,MINVAL_LEFT_RIGHT,outval);
blkout=outval;

%----------------------------------------------
%----------------------------------------------
function NC=print_NC(v,v_);
v=ShapeAsRow(v);
v_=v_(:);
% OLD !!!  NC=(v*v_) / (sqrt(sum(v.^2))*sqrt(sum(v_.^2)));
NC=(v*v_) /  sqrt( sum(v.^2) * sum(v_.^2) );
fprintf('Vector NC=%f\n',NC);
%v=ShapeAsRow(v);
%v_=ShapeAsRow(v_);
%meerNC=sum(v.*v_) / sqrt(sum(v.*v)*sum(v_.*v_));
%fprintf('meerNC=%f\n',meerNC);




%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lenaDashEC,insertedStream,insertedPayload,bpp9]=INSERT21(lena,th,PL);
global blkcountINSRT wmCount;
% notes:
% h=-2, 253 255 and 0 2, changeable, rest expandable! test with wm=0 and 1
% 255 255 pair, not expandable or changeable.

% classify
expandablePairMatrix=blkproc(lena,[2 1],@classifyExpandableFun,th);
changeablePairMatrix=blkproc(lena,[2 1],@classifyChangeableFun);
expandablePairArray=expandablePairMatrix';
expandablePairArray=ShapeAsRow(expandablePairArray(:));
changeablePairMatrix=changeablePairMatrix & ~expandablePairMatrix;
changeablePairArray=changeablePairMatrix';
changeablePairArray=ShapeAsRow(changeablePairArray(:));

totalNoOfCEpairs = sum(changeablePairArray) + sum(expandablePairArray);

hValueMatrix=blkproc(lena,[2 1],@getThe_h_DifferenceValuesFun);
Carray=getCarray(hValueMatrix,changeablePairMatrix); % store the original LSBs of changed pairs

imwrite(uint8(expandablePairArray*255),'exLM.pbm');
!pbmtojbg exLM.pbm exLM.jbg
fp=fopen('exLM.jbg','rb');
arrSymbols=fread(fp,'uchar');
fprintf('\nlength of location map=%d BYTES, [%d*8]+8=%d BITS [%d*9]+9=%d BITS\n', ...
	length(arrSymbols),length(arrSymbols),(length(arrSymbols)*8)+8,length(arrSymbols),(length(arrSymbols)*9)+9);
fclose(fp);
len9=[];
for i=1:length(arrSymbols),
	theValArr=dec2bin(arrSymbols(i),9);
	nineBin=double(uint8(theValArr)-48);
	len9=[len9 nineBin];
end
EOS=[1 1 1 1 1 1 1 1 1];
LarrComp=[len9 EOS]';  % locmap compressed as binary stream



% arrSymbols is compressed location map (theoretical just now!)
LC8 = length(arrSymbols)*8 + 8 + length(Carray); % amount of pairs needed for compressed LOCATION MAP and Carray (8bits)
LC9 = length(arrSymbols)*9 + 9 + length(Carray); % amount of pairs needed for compressed LOCATION MAP and Carray (9bits)
plBITS8 = totalNoOfCEpairs - LC8; % amount of pairs for embedding payload (no EOS, 8 bits)
plBITS9 = totalNoOfCEpairs - LC9; % amount of pairs for embedding payload (with EOS, 9 bits)
bpp8=plBITS8/prod(size(lena));
bpp9=plBITS9/prod(size(lena));
fprintf('\nbpp8=%.2f bpp9=%.2f\n',bpp8,bpp9);
fprintf('\nnoOfExpandablePairs=%d, numberOfChangeablePairs=%d\n',sum(expandablePairArray),sum(changeablePairArray));
fprintf('lengthPayload [8bits, no EOS]=%d, lengthPayload [9bits, with EOS]=%d\n',plBITS8,plBITS9);
fprintf('totalNoOfPairs=%d, TotalNoOfExpandableAndChangeablePairs=%d CarrayLength=%d\n',prod(size(lena))/2,totalNoOfCEpairs,length(Carray));
fprintf('length compressed location map = %d bits\n',length(LarrComp));
locmapsize=prod(size(lena))/2;

% concatenate the watermark ... LCP
wm=[LarrComp;Carray;PL(:)];

% insert changeable and expandable pairs
lenaDashEC=blkproc(lena,[2 1],@insertChangeableAndExpandable,changeablePairArray,expandablePairArray,wm);

%fprintf('\n\n\n\n\nmax lena=%.2f  min lena=%.2f\n',max(lena(:)),min(lena(:)));
%fprintf('max lenaDashEC=%.2f  min lenaDashEC=%.2f\n\n\n\n\n',max(lenaDashEC(:)),min(lenaDashEC(:)));

insertedStream = wm(1:totalNoOfCEpairs); % stream inserted into LSBs of changeable and expandable
%insertedPayload=insertedStream(length(Carray)+1:length(insertedStream));
lengthPL=totalNoOfCEpairs-length(LarrComp)-length(Carray);
insertedPayload=wm(length(LarrComp)+length(Carray)+1:totalNoOfCEpairs);

figure, GrayImage(expandablePairMatrix*255), title('Location map');



%%%%%%%%%%%%%%%%%%%% RECOVER %%%%%%%%%%%%%%%%%%%
% classify
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
function [LR,theBstreamRec,recoveredPayload]=RECOVER21(lenaDashEC);
global BstreamBlockCount BSTREAM  RIFCAEP_CNT CarrayCNT;
changeablePairMatrix2=blkproc(lenaDashEC,[2 1],@RECOVERclassifyChangeableFun);
changeablePairArray2=changeablePairMatrix2';
changeablePairArray2=ShapeAsRow(changeablePairArray2(:));

% FIND ALL CHANGEABLE PAIRS, store as B stream (BSTREAM)
dummy=blkproc(lenaDashEC,[2 1],@B_STREAM_recoverEC_LSBs,changeablePairArray2);
theBstreamRec=BSTREAM;

fprintf('length theBstreamRec=%d\n',length(theBstreamRec));

% expandable pairs are changeable!!
% would need to get the compressed locmap from B stream, decompress it, then 
% i have a map of where all the expanded pairs are. Here, it is just copied,
% assumed i have decompressed it
%LOCMAP=expandablePairMatrix;
%LOCMAP=LOCMAP';
%LOCMAP=ShapeAsRow(LOCMAP(:));
% see as.m and as2.m
asr=[];
thetest=0;
i=1;
%for i=1:9:length(len9),
while thetest == 0,
	get9=BSTREAM(i:i+8);
	kk=strread(num2str(get9),'%s')';
	kk=char(kk)';
	decval=double(bin2dec(kk));
	%fprintf('i=%d  decval=%d\n',i,decval);
	if decval == 511
		fprintf('EOS found\n');
		BSTREAM=BSTREAM(i+9:length(BSTREAM)); % chop BSTREAM (remove L, CP left)
		thetest=1;
		break;
	end
	asr=[asr decval];
	i=i+9;
end
fp2=fopen('lmrec.jbg','wb');
fwrite(fp2,asr,'uchar');
fclose(fp2);
!jbgtopbm lmrec.jbg lmrec.pbm
theLM=imread('lmrec.pbm');
LOCMAP=double(theLM);
fprintf('length BSTREAM=%d\n',length(BSTREAM));
fprintf('length LOCMAP=%d\n',length(LOCMAP));

% from all the changeable pairs, pick out the exapandable pairs
% changeablePairArray2 contains locations of all changeable and
% expandable pairs (expandable pairs are changeable!)
recoverExpandedArray = changeablePairArray2.*LOCMAP;
% pick out the pairs that were changed
recoverChangedArray=changeablePairArray2 .* ~LOCMAP;

% recover the original pair values that were changed and the original pair values that were expanded.
LR=blkproc(lenaDashEC,[2 1],@recoverImageFromChangeableAndExpandablePairs,recoverChangedArray,recoverExpandedArray,BSTREAM);
% stats
% CarrayCNT points to one place beyond length of Carray (see recoverImageFromChangeableAndExpandablePairs function).
%recoveredPayload=theBstreamRec(CarrayCNT:length(theBstreamRec)); 
recoveredPayload=BSTREAM(CarrayCNT:length(BSTREAM)); 

