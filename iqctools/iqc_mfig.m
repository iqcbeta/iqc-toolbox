function iqc_mfig(name1,name2)
% function iqc_mfig(name1,name2)
%
% creates/edits a drawing (segments, arrows, circles, texts, boxes etc.)
% the data is loaded from iqc_mfig_name1.mat and 
% saved in the form of iqc_mfig_name2.mat (ascii mat-files)
% outputs a LaTeX2e drawing script
%
% name1: name of the file (name1_iqc_mfig.mat) to be edited
%    default: test
% name2: the name under which the drawing should be saved
%    default: name1
%
% To SAVE AND EXIT: press "e"
% To SAVE: press "s"
% To DELETE A GRAPHICS OBJECT: press "d"; for each of the graphics object present
%                              a green dot will appear; click on the green dot to be
%                              deleted; sometimes two different objects will have
%                              coinciding green dots; then the younger object is
%                              eliminated first
% To ENTER A REFERENCE POINT: click left mouse button 
%                             [a small blue square indicates the position]
%     !!! ATTENTION !!! IQC_MFIG rounds any coordinate to an integer
%                       note the size of the IQC_MFIG grid is 2, not 1
% To ENTER A LINE: make sure a reference point is set; this will be one
%                  end of the line; click center mouse button
%                  where you want the other end of the line to be
%                  [line will appear on the screen, and its end will become
%                   the new reference point]
%     !!! ATTENTION !!! Any slope is automatically approximated by a LaTeX-supported one
% To ENTER A VECTOR: same as line, but use right mouse button
% To ENTER A TEXT: make sure a reference point is set; this will be the
%                  lower left corner of the text; press "t", then enter
%                  the text in the dialog box
% To ENTER A CIRCLE: the reference point will be the center, press "c",
%                    then click at a point on the circle; minimal radius is 1
% To ENTER A DISK: the reference point will be the center, press "k", 
%                  ("k" because a 2D disk is "kroog" in Russian)
%                   then click at a point on the circle; minimal radius is 1
% To ENTER A BOX (RECTANGULAR): the reference point will be one corner, 
%                 press "b", then click at the corner across the diagonal;
%                  you will be prompted for a text to be centered in the box              
% To ENTER A DASHBOX: same as for a box, but press "-" instead of "b"
% To ENTER AN OVAL: same as for a box, but press "o" instead of "b"
% To ENTER A BEZIER CURVE: the reference point will be one end of the curve;
%                 press "q" ("q" because of the "qbezier" command in LaTeX);
%                 then click at the middle reference point, and, finally, at the
%                 end point of the curve
% To ENTER A DOTTED BEZIER CURVE: same as for a Bezier curve, but using "."
%                 instead of "q"
% To ENTER A ROMB: the reference point will be one of the ends of an UPPER edge
%                 of the romb; press "r"; click on the other end of the edge
% To ENTER A DOTTED LINE: the reference point is one end of the line; press ",";
%                 click on the other end of the line


N=100; % horizontal size


if nargin<1,
   name1='test';
end
if nargin<2,
   name2=name1;   
end


% initialization
hf=figure;   % the drawing window
ha=axes('visible','off','XLim',[0 N],'YLim',[0 N]);

for k=0:2:N,   % horizontal grid
   line('XData',[0 N],'YData',[k k],'visible','on','Color','y');
end
for k=0:2:N,   % vertical grid
   line('XData',[k k],'YData',[0 N],'visible','on','Color','y');
end

% load the figure, if any
if eval(['exist(''iqc_mfig_' name1 '.mat'')'])==2,
   % a string "txt" and matrix "data" of size
   % k-by-m where k currently equals 8. 
   % A row of "data" corresponds to a graphics object
   % data(n,:)=[type x1 y1 x2 y2 t1 t2 z]
   % type = output of "ginput" when the command instruction is entered
   %     (for example, 2 for line, 3 for arrow, abs('t') for text, etc.)
   % (x1,y1) and (x2,y2) are always the "main coordinates" of the object
   % (t1,t2) -coordinates in "txt" of the text content of the object
   % z - extra parameter
   eval(['load iqc_mfig_' name1]);   
   % array to keep handles of graphics objects
   [k,m]=size(data);
   han=zeros(k,2); % to keep handles of graphics objects han(:,1)-lines, han(:,2)-text
   % display the picture
   for n=1:k,
      out=mfig_conv('graph',data(n,:),txt);
      gx=out{1}(1,:);
      gy=out{1}(2,:);
      txy=out{2};
      txto=out{3};
      stl=out{4};
      h1=line('XData',gx,'YData',gy,'LineStyle',stl,'visible','on','Color','r','MarkerSize',20);   
      h2=text('position',txy,'visible','on','string',txto);
      han(n,:)=[h1 h2];        
   end
else
   m=8;
   data=zeros(0,m); 
   txt='   ';
   k=0;
   han=zeros(0,2);
end


state=0;  % state of the system: 0 if no coordinate entered, 1 otherwise

 helpmes{1}='Mouse buttons: LB (left), CB (center), RB (right)';
 helpmes{2}='-------------------------------------------------';
 helpmes{3}='LINE:                LB  CB                      ';
 helpmes{4}='ARROW:               LB  RB                      ';
 helpmes{5}='TEXT:                LB  t                       ';
 helpmes{6}='CIRCLE:              LB  c  LB                   ';
 helpmes{7}='DISK:                LB  k  LB                   ';
 helpmes{8}='RECTANGULAR BOX:     LB  b  LB                   ';
 helpmes{9}='DASHED LINE BOX:     LB  -  LB                   ';
helpmes{10}='OVAL:                LB  o  LB                   ';
helpmes{11}='BEZIER CURVE:        LB  q  LB LB                ';
helpmes{12}='DOTTED BEZIER CURVE: LB  .  LB LB                ';
helpmes{13}='ROMB:                LB  r  LB                   ';
helpmes{14}='DOTTED LINE:         LB  ,  LB                   ';
helpmes{15}='-------------------------------------------------';
helpmes{16}='SAVE AND EXIT:       e                           ';
helpmes{17}='SAVE:                s                           ';
helpmes{18}='DELETE:              d   LB                      ';
helpmes{19}='HELP:                h                           ';

% drawing phase
while 1,
   [x,y,b]=ginput(1); % x,y-coordinates, b-key
   switch b,
   case 1, % enter the main coordinate
      % delete old marker, if any      
      if state==1,
          delete(hh);
      end
      state=1;
      x1=round(x);
      y1=round(y);
      % put marker
      xd=x1+[-0.2 -0.2 0.2 0.2 -0.2];
      yd=y1+[-0.2 0.2 0.2 -0.2 -0.2];
      hh=line('XData',xd,'YData',yd,'visible','on','Color','b','MarkerSize',20);
   case abs('h'),
      for khelp=1:length(helpmes),
         disp(helpmes{khelp})
      end
   case abs('d'),
      % delete old marker, if any      
      if state==1,
         delete(hh);
         state=0;
      end   
      hd=zeros(1,k);  % handles for "delete" marks
      for n=1:k,
         xd=data(n,2)+0.5*cos(0:(9*pi/10):100);
         yd=data(n,3)+0.5*sin(0:(9*pi/10):100);
         hd(n)=line('XData',xd,'YData',yd,'visible','on','Color','g','MarkerSize',20);
      end
      [x,y,z]=ginput(1);
      xd1=round(x);
      yd1=round(y);
      for n=k:-1:1,
         if (data(n,2)==xd1)&(data(n,3)==yd1),
            delete(han(n,1));
            delete(han(n,2));
            if (data(n,1)~=abs('q'))&(data(n,1)~=abs('.')),
               txt=[txt(1:data(n,6)-1) txt(data(n,7)+1:length(txt))];
               ldl=data(n,7)-data(n,6)+1; % length of the text deleted
            else
               ldl=0;
            end
            for n1=(n+1):k,
               if data(n1,1)~=abs('q'),
                  data(n1,6:7)=data(n1,6:7)-ldl;
               end
            end
            data=[data(1:n-1,:);data(n+1:k,:)];
            han=[han(1:n-1,:);han(n+1:k,:)];
            break
         end
      end   
      for n=1:k,
         delete(hd(n));
      end
      k=size(data,1);
   case abs('u'),
      if state==1,
         state=0;
         delete(hh);
      else
         if k>0,
            delete(han(k,1));
            delete(nah(k,2));
            txt=[txt(1:data(k,6)-1) txt(data(k,7)+1:length(txt))];
            k=k-1;
            han=han(1:k,:);
            data=data(1:k,:);
         end
      end
   case abs('e'),
      if state==1,
         delete(hh);
      end
      break;
   case abs('s'),
      eval(['save(''iqc_mfig_' name2 ''',''data'',''txt'')'])
   otherwise
      if state==1,
         out=mfig_conv('input',[b x1 y1 round(x) round(y)],txt);
         data=[data;out{1}];
         k=size(data,1);
         txt=[txt out{3}];
         delete(hh);         
         if ~isempty(out{2}),
            x1=out{2}(1);
            y1=out{2}(2);
            % put marker
            xd=x1+[-0.2 -0.2 0.2 0.2 -0.2];
            yd=y1+[-0.2 0.2 0.2 -0.2 -0.2];
            hh=line('XData',xd,'YData',yd,'visible','on','Color','b','MarkerSize',20);            
         else
            state=0;
         end
         if ~isempty(out{1}),   % display new object, if exists
            out=mfig_conv('graph',data(k,:),txt);
            gx=out{1}(1,:);
            gy=out{1}(2,:);
            txy=out{2};
            txto=out{3};
            stl=out{4};
            h1=line('XData',gx,'YData',gy,'LineStyle',stl,'visible','on','Color','r','MarkerSize',20);   
            h2=text('position',txy,'visible','on','string',txto, ...
               'HorizontalAlignment','center');
            han=[han;h1 h2];
         end
      end
   end
end

% saving txt and data into name_mfig.mat
eval(['save(''iqc_mfig_' name2 ''',''data'',''txt'')'])

latex=['\begin{figure}[h]'];
disp(latex)
latex=['\setlength{\unitlength}{0.09cm}'];
disp(latex)
% find box dimensions
xmin=N;
ymin=N;
xmax=0;
ymax=0;
for n=1:k,
   xmin=min([xmin data(n,2) data(n,4)]);
   xmax=max([xmax data(n,2) data(n,4)]);
   ymin=min([ymin data(n,3) data(n,5)]);
   ymax=max([ymax data(n,3) data(n,5)]);   
end
s1=num2str(xmax-xmin);
s2=num2str(ymax-ymin);
s3=num2str(xmin);
s4=num2str(ymin);
latex=['\begin{center}\begin{picture}(' s1 ',' s2 ')(' s3 ',' s4 ')'];
disp(latex)

for n=1:k,
   out=mfig_conv('latex',data(n,:),txt);
   for n1=1:length(out),
      disp(out{n1})
   end
end

latex='\end{picture}\end{center}';
disp(latex)
latex=['\caption{\ \label{fig:' name2 '}}'];
disp(latex)
latex='\end{figure}';
disp(latex)

close(hf)



function out=mfig_conv(ask,dat,txt)
% out=mfig_conv(ask,dat,txt)
%
% this function determines how iqc_mfig handles graphics objects
% 
% input argument "ask" determines the type of task:
%
% ask='input': input data for a graphics object 
%        out = {dato, refo, txto}
% ask='graph': graphics data for a graphics object described by a row vector dat
%        out = {gxy,txy,txto,stl}
% ask='latex': latex data for a graphics object described by a row vector dat
%        out = {line1,line2,...}
% dat = [type_of_entry, x1, y1, x2, y2, t1, t2, pos]
%       where type_of_entry = 2 for lines, 3 for arrows, abs(letter) for other types
%             x1,y1 - basic coordinates of the object
%             x2,y2 - secondary coordinates of the object
%             t1,t2 - begin and end of the text (no text if t2<t1)
%             pos   - a number from {1,2,...,9}, used for label positioning in boxes
%       in the 'input' mode, some of the components of "dat" are irrelevant

ntxt=length(txt);
stl='-';

switch dat(1),
case 2,
   switch ask,
   case 'input',
      [p,q,l]=lineslope(dat(4)-dat(2),dat(5)-dat(3));
      if l>0,
         dato=zeros(1,8);
         dato(1:3)=dat(1:3);
         if p~=0,
            dato(4:5)=[dat(2)+sign(p)*l dat(3)+q*l/abs(p)];
         else
            dato(4:5)=[dat(2) dat(3)+l*sign(q)];
         end
         dato(6:8)=[ntxt+1 ntxt 0];
         refo=dato(4:5);
         txto='';
      else
         dato=zeros(0,8);
         refo=dat(2:3);
         txto='';
      end  
      out={dato,refo,txto};
   case 'graph',         
      gxy=[dat(2) dat(4);dat(3) dat(5)];
      txy=[0 0];
      txto='';
      out={gxy,txy,txto,stl};
   case 'latex',
      s1=num2str(dat(2));
      s2=num2str(dat(3));
      [p,q,l]=lineslope(dat(4)-dat(2),dat(5)-dat(3));
      s3=num2str(p);
      s4=num2str(q);
      s5=num2str(l);
      out={['\put(' s1 ',' s2 '){\line(' s3 ',' s4 '){' s5 '}}']};        
   otherwise
      error('Inadmissible argument in mfig_conv')
   end   
case 3,
   switch ask,
   case 'input',
      [p,q,l]=vectorslope(dat(4)-dat(2),dat(5)-dat(3));
      if l>0,
         dato=zeros(1,8);
         dato(1:3)=dat(1:3);
         if p~=0,
            dato(4:5)=[dat(2)+sign(p)*l dat(3)+q*l/abs(p)];
         else
            dato(4:5)=[dat(2) dat(3)+l*sign(q)];
         end
         dato(6:8)=[ntxt+1 ntxt 0];
         refo=dato(4:5);
         txto='';
      else
         dato=zeros(0,8);
         refo=dat(2:3);
         txto='';
      end  
      out={dato,refo,txto};
   case 'graph',
      v=[dat(2)-dat(4);dat(3)-dat(5)];
      t=pi/10;
      v1=[cos(t) sin(t);-sin(t) cos(t)]*v/norm(v);
      v2=[cos(t) -sin(t);sin(t) cos(t)]*v/norm(v);         
      gxy=[dat(4)+[dat(2)-dat(4) 0 v1(1) v2(1) 0];dat(5)+[dat(3)-dat(5) 0 v1(2) v2(2) 0]];
      txy=[0 0];
      txto='';
      out={gxy,txy,txto,stl};
   case 'latex',
      s1=num2str(dat(2));
      s2=num2str(dat(3));
      [p,q,l]=vectorslope(dat(4)-dat(2),dat(5)-dat(3));
      s3=num2str(p);
      s4=num2str(q);
      s5=num2str(l);
      out={['\put(' s1 ',' s2 '){\vector(' s3 ',' s4 '){' s5 '}}']};        
   otherwise
      error('Inadmissible argument in mfig_conv')
   end
case abs('t'),
   switch ask,
   case 'input',
      z=inputdlg({'Text:'},'IQC_MFIG TEXT INPUT',1,{''});
      if isempty(z), z={''}; end
      txto=z{1};
      nz=length(txto);
      dato=[dat(1:5) ntxt+1 ntxt+nz 0];
      refo=zeros(0,2);
      out={dato,refo,txto};
   case 'graph',
      gxy=zeros(2,0);
      txy=[dat(2) dat(3)];
      txto=txt(dat(6):dat(7));
      out={gxy,txy,txto,stl};
   case 'latex',
      s1=num2str(dat(2));
      s2=num2str(dat(3));   
      s3=txt(dat(6):dat(7));
      out={['\put(' s1 ',' s2 '){$' s3 '$}']};
   otherwise
      error('Inadmissible argument in mfig_conv')
   end
case abs('c'),
   switch ask,
   case 'input',
      [x,y,z]=ginput(1);
      rad=round(sqrt((x-dat(2))^2+(y-dat(3))^2));
      if rad==0,
         rad=1;
      end
      dato=[abs('c') dat(2)-rad dat(3)-rad dat(2)+rad dat(3)+rad ntxt+1 ntxt 0];
      refo=zeros(0,2);
      txto='';
      out={dato,refo,txto};
   case 'graph',
      rad=(dat(4)-dat(2))/2;
      x1=(dat(4)+dat(2))/2;
      y1=(dat(5)+dat(3))/2;
      gxy=[x1+rad*cos(linspace(0,2*pi,100));y1+rad*sin(linspace(0,2*pi,100))];
      txy=[0 0];
      txto='';
      out={gxy,txy,txto,stl};
   case 'latex',
      s1=num2str((dat(2)+dat(4))/2);
      s2=num2str((dat(3)+dat(5))/2);
      s3=num2str(dat(4)-dat(2));
      out={['\put(' s1 ',' s2 '){\circle{' s3 '}}']};       
   otherwise
      error('Inadmissible argument in mfig_conv')
   end   
case abs('k'),
   switch ask,
   case 'input',
      [x,y,z]=ginput(1);
      rad=round(sqrt((x-dat(2))^2+(y-dat(3))^2));
      if rad==0,
         rad=1;
      end
      dato=[abs('k') dat(2)-rad dat(3)-rad dat(2)+rad dat(3)+rad ntxt+1 ntxt 0];
      refo=zeros(0,2);
      txto='';
      out={dato,refo,txto};
   case 'graph',
      rad=(dat(4)-dat(2))/2;
      x1=(dat(4)+dat(2))/2;
      y1=(dat(5)+dat(3))/2;
      gxy=[x1+rad*cos(linspace(0,200*pi,211));y1+rad*sin(linspace(0,200*pi,211))];
      txy=[0 0];
      txto='';
      out={gxy,txy,txto,stl};
   case 'latex',
      s1=num2str((dat(2)+dat(4))/2);
      s2=num2str((dat(3)+dat(5))/2);
      s3=num2str(dat(4)-dat(2));
      out={['\put(' s1 ',' s2 '){\circle*{' s3 '}}']};       
   otherwise
      error('Inadmissible argument in mfig_conv')
   end  
case abs('b'),
   switch ask,
   case 'input',
      [x,y,z]=ginput(1);
      x=round(x);
      y=round(y);
      x1=min(x,dat(2));
      y1=min(y,dat(3));
      x2=max(x,dat(2));
      y2=max(y,dat(3));
      if x1==x2,
         x1=x1-1;
         x2=x2+1;
      end
      if y1==y2,
         y1=y1-1;
         y2=y2+1;
      end
      z=inputdlg({'Text:'},'IQC_MFIG BOX TEXT INPUT',1,{''});
      if isempty(z), z={''}; end
      txto=z{1};
      nz=length(txto);     
      dato=[abs('b') x1 y1 x2 y2 ntxt+1 ntxt+nz 0];
      refo=zeros(0,2);
      out={dato,refo,txto};
   case 'graph',
      gxy=[dat(2) dat(2) dat(4) dat(4) dat(2);dat(3) dat(5) dat(5) dat(3) dat(3)];
      txy=[0.5*(dat(2)+dat(4)) (dat(3)+dat(5))/2];
      txto=txt(dat(6):dat(7));
      out={gxy,txy,txto,stl};
   case 'latex',
      s1=num2str(dat(2));
      s2=num2str(dat(3));
      s3=num2str(dat(4)-dat(2));
      s4=num2str(dat(5)-dat(3));
      s5=txt(dat(6):dat(7));
      out={['\put(' s1 ',' s2 '){\framebox(' s3 ',' s4 '){$ ' s5 '$}}']};       
   otherwise
      error('Inadmissible argument in mfig_conv')
   end    
case abs('-'),
   switch ask,
   case 'input',
      [x,y,z]=ginput(1);
      x=round(x);
      y=round(y);
      x1=min(x,dat(2));
      y1=min(y,dat(3));
      x2=max(x,dat(2));
      y2=max(y,dat(3));
      if x1==x2,
         x1=x1-1;
         x2=x2+1;
      end
      if y1==y2,
         y1=y1-1;
         y2=y2+1;
      end
      z=inputdlg({'Text:'},'IQC_MFIG DASHBOX TEXT INPUT',1,{''});
      if isempty(z), z={''}; end
      txto=z{1};
      nz=length(txto);     
      dato=[abs('-') x1 y1 x2 y2 ntxt+1 ntxt+nz 0];
      refo=zeros(0,2);
      out={dato,refo,txto};
   case 'graph',
      gxy=[dat(2) dat(2) dat(4) dat(4) dat(2);dat(3) dat(5) dat(5) dat(3) dat(3)];
      txy=[0.5*(dat(2)+dat(4)) (dat(3)+dat(5))/2];
      txto=txt(dat(6):dat(7));
      out={gxy,txy,txto,'--'};
   case 'latex',
      s1=num2str(dat(2));
      s2=num2str(dat(3));
      s3=num2str(dat(4)-dat(2));
      s4=num2str(dat(5)-dat(3));
      s5=txt(dat(6):dat(7));
      out={['\put(' s1 ',' s2 '){\dashbox(' s3 ',' s4 '){$ ' s5 '$}}']};       
   otherwise
      error('Inadmissible argument in mfig_conv')
   end    
case abs('o'),
   switch ask,
   case 'input',
      [x,y,z]=ginput(1);
      x=round(x);
      y=round(y);
      x1=min(x,dat(2));
      y1=min(y,dat(3));
      x2=max(x,dat(2));
      y2=max(y,dat(3));
      if x1==x2,
         x1=x1-1;
         x2=x2+1;
      end
      if y1==y2,
         y1=y1-1;
         y2=y2+1;
      end
      z=inputdlg({'Text:'},'IQC_MFIG OVAL TEXT INPUT',1,{''});
      if isempty(z), z={''}; end
      txto=z{1};
      nz=length(txto);     
      dato=[abs('o') x1 y1 x2 y2 ntxt+1 ntxt+nz 0];
      refo=zeros(0,2);
      out={dato,refo,txto};
   case 'graph',
      tt=linspace(0,pi/2,50);
      cost=cos(tt);
      sint=sin(tt);
      rad=0.4*min(dat(4)-dat(2),dat(5)-dat(3));
      gx=[dat(2) dat(2)+rad*(1-cost) dat(4)-rad*(1-sint) dat(4)-rad*(1-cost) dat(2)+rad*(1-sint)];
      gy=[dat(3)+rad dat(5)-rad*(1-sint) dat(5)-rad*(1-cost) dat(3)+rad*(1-sint) dat(3)+rad*(1-cost)];
      gxy=[gx;gy];
      txy=[0.5*(dat(2)+dat(4)) (dat(3)+dat(5))/2];
      txto=txt(dat(6):dat(7));
      out={gxy,txy,txto,stl};
   case 'latex',
      s1=num2str(0.5*(dat(2)+dat(4)));
      s2=num2str(0.5*(dat(3)+dat(5)));
      s3=num2str(dat(4)-dat(2));
      s4=num2str(dat(5)-dat(3));
      s5=txt(dat(6):dat(7));
      out{1}=['\put(' s1 ',' s2 '){\oval(' s3 ',' s4 ')}'];   
      out{2}=['\put(' s1 ',' s2 '){\makebox(0,0){' s5 '}}'];
   otherwise
      error('Inadmissible argument in mfig_conv')
   end  
case abs('q'),
   switch ask,
   case 'input',
      % get the second point
      [x,y,z]=ginput(1);
      x=round(x);
      y=round(y);
      % put marker for the second reference point
      xd=x+[-0.2 -0.2 0.2 0.2 -0.2 0 dat(2)-x];
      yd=y+[-0.2 0.2 0.2 -0.2 -0.2 0 dat(3)-y];
      hh1=line('XData',xd,'YData',yd,'visible','on','Color','b','MarkerSize',20);    
      % get the third point
      [x1,y1,z1]=ginput(1);
      delete(hh1);
      dato=[abs('q') dat(2:3) x y round(x1) round(y1) 0];
      refo=[round(x1) round(y1)];
      txto='';
      out={dato,refo,txto};
   case 'graph',
      tt=linspace(0,pi/2,200);
      cost=cos(tt);
      sint=sin(tt);
      v1=[dat(2)-dat(4);dat(3)-dat(5)];
      v2=[dat(6)-dat(4);dat(7)-dat(5)];
      v0=[dat(4);dat(5)];
      gxy=(v0+v1+v2)*ones(1,length(cost))-v1*cost-v2*sint;
      txy=[0 0];
      txto='';
      out={gxy,txy,txto,stl};
   case 'latex',
      s1=num2str(dat(2));
      s2=num2str(dat(3));
      s3=num2str(dat(4));
      s4=num2str(dat(5));
      s5=num2str(dat(6));
      s6=num2str(dat(7));
      out={['\qbezier(' s1 ',' s2 ')(' s3 ',' s4 ')(' s5 ',' s6 ')']};
   otherwise
      error('Inadmissible argument in mfig_conv')
   end  
case abs('.'),
   switch ask,
   case 'input',
      % get the second point
      [x,y,z]=ginput(1);
      x=round(x);
      y=round(y);
      % put marker for the second reference point
      xd=x+[-0.2 -0.2 0.2 0.2 -0.2 0 dat(2)-x];
      yd=y+[-0.2 0.2 0.2 -0.2 -0.2 0 dat(3)-y];
      hh1=line('XData',xd,'YData',yd,'visible','on','Color','b','MarkerSize',20);    
      % get the third point
      [x1,y1,z1]=ginput(1);
      delete(hh1);
      dato=[abs('.') dat(2:3) x y round(x1) round(y1) 0];
      refo=[round(x1) round(y1)];
      txto='';
      out={dato,refo,txto};
   case 'graph',
      tt=linspace(0,pi/2,200);
      cost=cos(tt);
      sint=sin(tt);
      v1=[dat(2)-dat(4);dat(3)-dat(5)];
      v2=[dat(6)-dat(4);dat(7)-dat(5)];
      v0=[dat(4);dat(5)];
      gxy=(v0+v1+v2)*ones(1,length(cost))-v1*cost-v2*sint;
      txy=[0 0];
      txto='';
      out={gxy,txy,txto,':'};
   case 'latex',
      s1=num2str(dat(2));
      s2=num2str(dat(3));
      s3=num2str(dat(4));
      s4=num2str(dat(5));
      s5=num2str(dat(6));
      s6=num2str(dat(7));
      out={['\qbezier[30](' s1 ',' s2 ')(' s3 ',' s4 ')(' s5 ',' s6 ')']};
   otherwise
      error('Inadmissible argument in mfig_conv')
   end  
case abs('r'),
   switch ask,
   case 'input',
      [x,y,z]=ginput(1);
      x=round(x);
      y=round(y);
      [p,q,l]=lineslope(x-dat(2),y-dat(3));
      if l*p*q==0,
         dato=zeros(0,8);
         refo=dat(2:3);
         txto='';
      else
         x=dat(2)+sign(p)*l;
         y=dat(3)+q*l/abs(p);
         if p*q>0,
            x1=min(x,dat(2));
            x2=x1+2*abs(x-dat(2));
            y2=max(y,dat(3));
            y1=y2-2*abs(y-dat(3));
         else
            x2=max(x,dat(2));
            x1=x2-2*abs(x-dat(2));
            y2=max(y,dat(3));
            y1=y2-2*abs(y-dat(3));            
         end   
         z=inputdlg({'Text:'},'IQC_MFIG ROMB TEXT INPUT',1,{''});
         if isempty(z), z={''}; end
      	txto=z{1};
      	nz=length(txto);     
      	dato=[abs('r') min(x1,x2) min(y1,y2) max(x1,x2) max(y1,y2) ntxt+1 ntxt+nz 0];
      	refo=zeros(0,2);
         out={dato,refo,txto};
      end
   case 'graph',
      x=(dat(2)+dat(4))/2;
      y=(dat(3)+dat(5))/2;
      gxy=[dat(2) x dat(4) x dat(2);y dat(5) y dat(3) y];
      txy=[0.5*(dat(2)+dat(4)) (dat(3)+dat(5))/2];
      txto=txt(dat(6):dat(7));
      out={gxy,txy,txto,stl};
   case 'latex',
      x=(dat(2)+dat(4))/2;
      y=(dat(3)+dat(5))/2;      
      s1=num2str(dat(2));
      s2=num2str(dat(3));
      s3=num2str(dat(4));
      s4=num2str(dat(5));
      s5=txt(dat(6):dat(7));
      sx=num2str(x);
      sy=num2str(y);
      [p,q,l]=lineslope(x-dat(2),y-dat(3));
      sp=num2str(p);
      sq=num2str(q);
      sl=num2str(l);
      out{1}=['\put(' s1 ',' sy '){\line(' sp ',' sq '){' sl '}}'];
      out{2}=['\put(' s3 ',' sy '){\line(-' sp ',' sq '){' sl '}}'];
      out{3}=['\put(' s1 ',' sy '){\line(' sp ',-' sq '){' sl '}}'];
      out{4}=['\put(' s3 ',' sy '){\line(-' sp ',-' sq '){' sl '}}'];
      out{5}=['\put(' sx ',' sy '){\makebox(0,0){' s5 '}}'];
   otherwise
      error('Inadmissible argument in mfig_conv')
   end   
case abs(','),
   switch ask,
   case 'input',
      % get the second point
      [x,y,z]=ginput(1);
      x=round(x);
      y=round(y);
      dato=[abs(',') dat(2:3) x y ntxt+1 ntxt 0];
      refo=[x y];
      txto='';
      out={dato,refo,txto};
   case 'graph',
      gxy=[dat(2) dat(4);dat(3) dat(5)];
      txy=[0 0];
      txto='';
      out={gxy,txy,txto,':'};
   case 'latex',
      s1=num2str(dat(2));
      s2=num2str(dat(3));
      s3=num2str(dat(4));
      s4=num2str(dat(5));
      s5=num2str((dat(2)+dat(4))/2);
      s6=num2str((dat(3)+dat(5))/2);
      out={['\qbezier[30](' s1 ',' s2 ')(' s5 ',' s6 ')(' s3 ',' s4 ')']};
   otherwise
      error('Inadmissible argument in mfig_conv')
   end   
end




function [po,qo,lo]=vectorslope(a,b)
% determines the LaTeX admissible vector slope which is the closest to (a,b)

p=zeros(1,24);
q=zeros(1,24);
p(1:13)=[1 0 1 1 1 1 2 2 3 3 3 4 4];
q(1:13)=[0 1 1 2 3 4 1 3 1 2 4 1 3];
p(14:24)=p(3:13);
q(14:24)=-q(3:13);

lpq=p.^2+q.^2;
s=(b*p-a*q)./lpq;
l=sign(a*p+b*q).*round(sqrt((a+s.*q).^2+(b-s.*p).^2)./sqrt(lpq));

[y,ind]=min((a-p.*l).^2+(b-q.*l).^2);

po=p(ind);
qo=q(ind);
lo=l(ind);
if po==0,
   lo=qo*lo;
else
   lo=po*lo;
end
if lo<0,
   lo=-lo;
   po=-po;
   qo=-qo;
end

function [po,qo,lo]=lineslope(a,b)
% determines the LaTeX admissible vector slope which is the closest to (a,b)

p=zeros(1,48);
q=zeros(1,48);
p(1:25)=[1 0   1 1 1 1 1 1 2 2 2 3   3 3 3 4 4 4 5 5 5 5   5 6 6];
q(1:25)=[0 1   1 2 3 4 5 6 1 3 5 1   2 4 5 1 3 5 1 2 3 4   6 1 5];
p(26:48)=p(3:25);
q(26:48)=-q(3:25);

lpq=p.^2+q.^2;
s=(b*p-a*q)./lpq;
l=sign(a*p+b*q).*round(sqrt((a+s.*q).^2+(b-s.*p).^2)./sqrt(lpq));

[y,ind]=min((a-p.*l).^2+(b-q.*l).^2);

po=p(ind);
qo=q(ind);
lo=l(ind);
if po==0,
   lo=qo*lo;
else
   lo=po*lo;
end
if lo<0,
   lo=-lo;
   po=-po;
   qo=-qo;
end



      
      
      
