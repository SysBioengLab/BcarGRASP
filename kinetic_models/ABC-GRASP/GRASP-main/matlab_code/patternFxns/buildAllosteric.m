function buildAllosteric(metList,reactionName,negEffectors,posEffectors)
% Builds allosteric reaction file.
%
%
% USAGE:
%
%    buildAllosteric(metList, reactionName, negEffectors, posEffectors)
%
% INPUT:
%    metList (struct):      list with all mets concentrations showing up in the pseudo-first-order rate constants 
%    reactionName (char):   reaction name
%    negEffectors (vector): list with negative effectors
%    posEffectors (vector): list with positive effectors
%
% OUTPUT:
%    written .m file with the reaction mechanism
%
% .. Authors:
%       - Pedro Saa     2016 original code 

% 1. Get output file handler
currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
filepath = fullfile(currentPath{1}, '..', '..', 'temp', 'reactions', [reactionName,'.m']);

fid = fopen(filepath, 'w'); 
if fid == -1
    error(['File not found: ', filepath, ...
            newline, ...
           'Please make sure the folder ', fullfile(currentPath{1}, '..', '..', 'temp', 'reactions'), ...
           ' exists.'])
end

% 2. Write initial parameters
c = '%';
if (isempty(metList))
    fprintf(fid,['function v = ',reactionName,'(X,negEff,posEff,Kr,KposEff,KnegEff,L,n) \n']);
end
if(~isempty(metList))
    if ~isempty(negEffectors) && ~isempty(posEffectors)
        fprintf(fid,['function v = ',reactionName,'(X,negEff,posEff,Kr,KnegEff,KposEff,L,n) \n']);   
    elseif ~isempty(negEffectors)
        fprintf(fid,['function v = ',reactionName,'(X,negEff,Kr,KnegEff,L,n) \n']);   
    elseif ~isempty(posEffectors)
        fprintf(fid,['function v = ',reactionName,'(X,posEff,Kr,KposEff,L,n) \n']);
    else
        fprintf(fid,['function v = ',reactionName,'(X,Kr,L,n) \n']);
    end        
end
fprintf(fid,'%s Parameters definition \n',c);
strform = strcat(reactionName,'Catalytic');

% Active form reaction rate
fprintf(fid,['[vR,eR] = ',strform,'(X,Kr); \n']);

% 3. Print reaction terms
fprintf(fid,'Q = L*eR.^n; \n');
if ~isempty(negEffectors)    
    fprintf(fid,'KnegEff = KnegEff(ones(size(negEff,2),1),:); \n');
    fprintf(fid,'Q = Q.*((1 + sum(negEff''./KnegEff,2)).^n)''; \n');
end
if ~isempty(posEffectors)
    fprintf(fid,'KposEff = KposEff(ones(size(posEff,2),1),:); \n');
    fprintf(fid,'Q = Q.*((1 + sum(posEff''./KposEff,2)).^-n)''; \n');
end
fprintf(fid,'%s Reaction rate \n',c);
fprintf(fid,'v = n*vR./(1 + Q);');
fclose(fid);