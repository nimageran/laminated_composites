function z=ModulFit(GA_thetaHalf)

GA_thetaHalf=round(GA_thetaHalf);
display(GA_thetaHalf);
nnn=numel(GA_thetaHalf); % Number of GA thetaHalf 
THETA=zeros(1,2*nnn);
for j=1:nnn
    if(GA_thetaHalf(j))==1;
        THETA(1,2*j-1:2*j)=[0 0];
    elseif(GA_thetaHalf(j))==2;
       THETA(1,2*j-1:2*j)=[+45 -45];
    else
       THETA(1,2*j-1:2*j)=[90 90];
    end
end
display(GA_thetaHalf);
GA_thetaHalf=THETA

global NFE;
if isempty(NFE)
    NFE=0;
end


% Updating NFE number
NFE=NFE+1;

%%%%%% Write variables into TEXT2.TXT file 
dlmwrite('TEXT2.txt',GA_thetaHalf,'\t')

%%%%%% Writing the outputs by running the python code
disp('Start python')
system(['abaqus cae nogui=pythonScript.py'])
disp('Python is finished')

newpath=['D:\optLaminatedComp\JOB\',mat2str(NFE)];

loadFactor=importdata(fullfile( newpath,'loadFactor.txt'));
z=-loadFactor;
disp(class(z))
disp(z)
end