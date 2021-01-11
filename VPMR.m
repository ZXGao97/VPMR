function varargout = VPMR(varargin)
% VPMR MATLAB code for VPMR.fig
%      VPMR, by itself, creates a new VPMR or raises the existing
%      singleton*.
%
%      H = VPMR returns the handle to a new VPMR or the handle to
%      the existing singleton*.
%
%      VPMR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VPMR.M with the given input arguments.
%
%      VPMR('Property','Value',...) creates a new VPMR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VPMR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VPMR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VPMR

% Last Modified by GUIDE v2.5 11-Jan-2021 20:46:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VPMR_OpeningFcn, ...
                   'gui_OutputFcn',  @VPMR_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before VPMR is made visible.
function VPMR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VPMR (see VARARGIN)

% Choose default command line output for VPMR
handles.output = hObject;
handles.popindex=1;
handles.err=[];
handles.XX=0;
handles.pp=[];
handles.ww=[];
handles.Anew=[];
handles.Bnew=[];
handles.Cnew=[];
handles.radioindex=0;
handles.radioindex2=0;
handles.digitss=300;
handles.xmax=0;
handles.xmin=0;
handles.n=0;
handles.stepsize=0;
handles.sclae=0;
handles.errors=0;
handles.errs=[];
handles.s=[];
handles.accindex=1;
handles.mrterms=0;
handles.intindex=2;
handles.records={'Record'};
set(handles.pushbutton1,'Visible','Off');
set(handles.pushbutton3,'Visible','Off');
set(handles.pushbutton4,'Visible','Off');
set(handles.pushbutton6,'Visible','Off');
set(handles.text3,'Visible','Off');
set(handles.text4,'Visible','Off');
set(handles.text5,'Visible','Off');
set(handles.edit1,'Visible','Off');
set(handles.edit2,'Visible','Off');
set(handles.edit3,'Visible','Off');
set(handles.text7,'Visible','Off');
set(handles.text8,'Visible','Off');
set(handles.text11,'Visible','Off');
set(handles.text13,'Visible','Off');
set(handles.edit5,'Visible','Off');
set(handles.text15,'Visible','Off');
set(handles.edit6,'Visible','Off');
set(handles.text16,'Visible','Off');
set(handles.text19,'Visible','Off');
set(handles.text20,'Visible','Off');
set(handles.edit7,'Visible','Off');
set(handles.edit8,'Visible','Off');
set(handles.edit9,'Visible','Off');
set(handles.edit10,'Visible','Off');
set(handles.edit11,'Visible','Off');
set(handles.text22,'Visible','Off');
set(handles.text23,'Visible','Off');
set(handles.text24,'Visible','Off');
set(handles.text25,'Visible','Off');
set(handles.text26,'Visible','Off');
set(handles.popupmenu4,'Visible','Off');
set(handles.popupmenu5,'Visible','Off');
set(handles.popupmenu6,'Visible','Off');
set(handles.text27,'Visible','Off');
set(handles.radiobutton2,'Visible','Off');
set(handles.edit13,'Visible','Off');
set(handles.edit14,'Visible','Off');
set(handles.edit15,'Visible','Off');
set(handles.text28,'Visible','Off');
set(handles.text29,'Visible','Off');
set(handles.text30,'Visible','Off');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VPMR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VPMR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
%SOG---MAIN
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.popindex=get(handles.popupmenu1,'Value');
if handles.popindex==2 || handles.popindex==6
    acc=handles.accindex;
    x=str2double(get(handles.edit3,'string'));handles.xmax=x;
    n=str2double(get(handles.edit1,'string'));handles.n=n;
    if handles.popindex==2
        v=eval(get(handles.edit8,'string'));handles.v=v;
        stepsize=eval(get(handles.edit9,'string'));handles.stepsize=stepsize;
    else
        v=0;stepsize=0.002;
        handles.xmin=v;handles.stepsize=stepsize;
    end
    scale=str2double(get(handles.edit2,'string'));handles.scale=scale;
    s=get(handles.edit4,'string');handles.s=s;handles.records=[handles.records,{s}];
    set(handles.popupmenu7,'String',handles.records);
    set(handles.text8,'string','Now Loading');
    if handles.radioindex==0
        handles.err=exptexting2(x,n,scale,s,v,stepsize,acc);
    else
        s1=s;s2=get(handles.edit10,'string');
        sp=eval(get(handles.edit11,'string'));
        handles.err=exptexting5(x,n,scale,s1,s2,sp,v,stepsize,acc);
    end
    digitts=max(handles.err(3,:));
    set(handles.text8,'string',num2str(digitts));
    axes(handles.axes1);
    plot(handles.err(1,:),handles.err(3,:));
    xlabel('x');ylabel('Error');
    guidata(hObject, handles); 
elseif handles.popindex==4 || handles.popindex==5 || handles.popindex==8 || handles.popindex==9
    u=eval(get(handles.edit3,'string'));handles.xmax=u;
    n=eval(get(handles.edit1,'string'));handles.n=n;
    zz=eval(get(handles.edit2,'string'));handles.scale=zz;
    s=get(handles.edit4,'string');handles.s=[s,'--SOG'];handles.records=[handles.records,{s}];
    set(handles.popupmenu7,'String',handles.records);
    g=@(x) eval(s);
    if handles.popindex==4 || handles.popindex==5
        digitts=eval(get(handles.edit7,'string'));
        handles.digitts=digitts;
    else
        digitts=300;
    end
    mp.Digits(digitts);
    scale=n*zz;
    set(handles.text12,'string','Computing Integrals');drawnow;
    if handles.radioindex==0
        A=zeros(1,2*n);
        if handles.intindex==2
            for i=0:2*n-1
                t=i;
                set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
                A(i+1)=quadgk(q,0,pi,'AbsTol',handles.accindex,'MaxIntervalCount',1e9,'RelTol',0);
            end
        elseif handles.intindex==3
            w=handles.accindex;
            ww=ceil(sqrt(w));
            d=2;
            S=2*ww/d;
            N=100;
            [b1,b2]=grule(N);
            for i=0:2*n-1
                t=i;
                set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
                for tt=1:S
                    A(i+1)=A(i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
                end
                A(i+1)=A(i+1)/ww;
            end
        end
    else
        s1=s;s2=get(handles.edit10,'string');
        sp=eval(get(handles.edit11,'string'));
        g1=@(x) eval(s1);g2=@(x) eval(s2);sps=acos(2*exp(-sp*sp/scale)-1);
        A=zeros(1,2*n);
        if handles.intindex==2
            for i=0:2*n-1
                t=i;
                set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                q1=@(x) (2/pi)*cos(t*x).*g1(sqrt(-scale*log((1+cos(x))/2)));
                q2=@(x) (2/pi)*cos(t*x).*g2(sqrt(-scale*log((1+cos(x))/2)));
                A(i+1)=quadgk(q1,0,sps,'AbsTol',handles.accindex,'MaxIntervalCount',1e9,'RelTol',0)+quadgk(q2,sps,pi,'AbsTol',1.0e-12,'MaxIntervalCount',1e9,'RelTol',0);
            end
        elseif handles.intindex==3
            w=handles.accindex;
            ww=ceil(sqrt(w));
            d=2;
            S=2*ww/d;
            N=100;
            [b1,b2]=grule(N);
            for i=0:2*n-1
                t=i;
                set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                q1=@(x) (2/pi)*cos(t*x).*g1(sqrt(-scale*log((1+cos(x))/2)));
                q2=@(x) (2/pi)*cos(t*x).*g2(sqrt(-scale*log((1+cos(x))/2)));
                for tt=1:S
                    A(i+1)=A(i+1)+q1(((2*tt-1+b1-ww)/ww)*sps/2+sps/2)*b2'*sps/2;
                    A(i+1)=A(i+1)+q2(((2*tt-1+b1-ww)/ww)*(pi-sps)/2+(sps+pi)/2)*b2'*(pi-sps)/2;
                end
                A(i+1)=A(i+1)/ww;
            end
        end
    end
    A=mp(A,digitts);
    maxC=0.00;
    set(handles.text12,'string','Computing Linear Coefficients(VP)');drawnow;
    for r=0:2*n-1
        set(handles.text11,'string',[num2str(r/2/n*100),'%(2/12)']);drawnow;
        if(r==0)
            sum1=mp('0.00');
            for k=1:n
                sum1=mp(sum1+mp(-1)^mp(k)*mp(A(k+1)),digitts);
            end
            sum2=mp('0.00');
            for k=1:n-1
                sum2=mp(sum2+mp(-1)^mp(n+k)*(1-mp(k)/n)*mp(A(n+k+1)),digitts);
            end
            C(r+1)=mp(A(r+1)/2+sum1+sum2,digitts);
        end
        if(1<=r&&r<=n)
            sum1=mp('0.00');
            sum2=mp('0.00');
            for k=r:n
                sum1=mp(sum1+mp((-1)^(k-r))*mp(2*mp(k)/mp(k+r))* mp(factorial(mp(k+r))) / mp(factorial(mp(2*r))*factorial(mp(k-r)))*mp((mp(2)^mp(2*r-1)))*mp(A(k+1)),digitts );
            end
            for k=1:n-1
                sum2=mp(sum2+mp((-1)^(n+k-r))*mp(1-mp(k)/n)*2*mp((n+k)*(mp(1)/(n+k+r)))*(mp(factorial(mp(n+k+r))))/((mp(factorial(mp(2*r))))*(mp(factorial(mp(n+k-r)))))*mp(2^mp(2*r-1))*mp(A(n+k+1)),digitts);
            end
            C(r+1)=mp((sum1)+(sum2),digitts);
        end
        if(r>=n+1)
            sum1=mp('0.00');
            for k=(r-n:n-1)
                sum1=mp(sum1+mp(mp(-1)^mp(n+k-r))*mp(1-mp(k)/n)*2*mp(mp(n+k)/(mp(n+k+r)))*(mp(factorial(mp(n+k+r))))/(mp(factorial(mp(2*r)))*mp(factorial(mp(n+k-r))))*(mp(2^mp(2*r-1)))*mp(A(n+k+1)),digitts);
            end
            C(r+1)=mp(sum1,digitts);
        end
    end
    if handles.popindex==4 || handles.popindex==5
        v=eval(get(handles.edit8,'string'));handles.xmin=v;
        stepsize=eval(get(handles.edit9,'string'));handles.stepsize=stepsize;
        x=(v+stepsize):stepsize:u;
    else
        stepsize=0.002;v=0;
        handles.xmin=v;handles.stepsize=stepsize;
        if u>=10
            stepsize=0.2;
        end
        x=stepsize:stepsize:u;
    end
    if handles.radioindex2==0
        NNN=size(x,2);
%         for i=1:NNN
%             i;
%             a(i)=mp('0.00');
%         end
        set(handles.text12,'string','Testing the VP sum');drawnow;
        set(handles.text11,'string','(3/12)');drawnow;
%         for k=0:2*n-1
%             set(handles.text11,'string',[num2str(k/2/n*100),'%(3/12)']);drawnow;
%             for i=1:NNN
%                 y(i)=mp(exp(mp(-mp(k)*(mp(x(i))*mp(x(i)))/mp(scale))));
%                 a(i)=mp(a(i)+mp(C(k+1))*mp(y(i)),digitts);
%             end 
%         end
        CC=C.';xx=mp(x).^2;KK=mp(0:(2*n-1))/mp(scale);
        a=sum(CC.*exp(-KK'.*xx),1);
        eta=double(max(mp(abs(mp(a-(g(x)))))));
        if handles.radioindex==1
            eta1=[ones(1,round((sp-v)/stepsize)) zeros(1,round(NNN-(sp-v)/stepsize))];
            eta2=[zeros(1,round((sp-v)/stepsize)) ones(1,round(NNN-(sp-v)/stepsize))];
            gg=g1(x).*eta1+g2(x).*eta2;
            eta=double(max(mp(abs(mp(a-gg)))));
        end
        set(handles.text8,'string',[num2str(eta) '(VP)']);drawnow;
    end
    X=C;handles.XX=double(X(1));
    %%Model Reduction
    set(handles.text12,'string','Model Reduction');drawnow;
    set(handles.text11,'string','(4/12)');drawnow;
    n=2*n;mp.Digits(digitts);
    if handles.popindex==4 || handles.popindex==5
        v=eval(get(handles.edit8,'string'));
        stepsize=eval(get(handles.edit9,'string'));
        x=(v+stepsize):stepsize:u;
    else
        stepsize=0.002;
        if u>=10
            stepsize=0.2;
        end
        x=stepsize:stepsize:u;
    end
    N=size(x,2);
    y=g(x);
    if handles.radioindex==1
        NNN=size(x,2);
        eta1=[ones(1,round((sp-v)/stepsize)) zeros(1,round(NNN-(sp-v)/stepsize))];
        eta2=[zeros(1,round((sp-v)/stepsize)) ones(1,round(NNN-(sp-v)/stepsize))];
        y=g1(x).*eta1+g2(x).*eta2;
    end
    y=y.';
    A=mp(-mp(diag(1:n-1))/mp(n*zz/2),digitts);
    A2=mp(-mp((1:n-1))/mp(n*zz/2),digitts);
    B=mp(zeros(n-1,1),digitts);
    C=mp(zeros(1,n-1),digitts);
    for i=1:n-1
        if(X(i+1)>0)
            C(i)=mp(sqrt(mp(X(i+1),digitts)),digitts);
            B(i)=mp(sqrt(mp(X(i+1),digitts)),digitts);
        else
            C(i)=mp(sqrt(mp(-X(i+1),digitts)),digitts);
            B(i)=mp(-sqrt(mp(-X(i+1),digitts)),digitts);
        end
    end
    BBB=B*B.';
    CCC=C.'*C;
    AA=A2+A2.';
    set(handles.text12,'string','Preparing P and Q');drawnow;
    set(handles.text11,'string',' ');drawnow;
    P=-BBB./AA;Q=-CCC./AA;
    set(handles.text12,'string','Cholesky');drawnow;
    mp.Digits(digitts);
    for i=1:n-1
        set(handles.text11,'string',[num2str(i/n*100),'%(6/12)']);drawnow;
        P(i,i)=sqrt(P(i,i)-sum(P(i,1:i-1).^2));
        Q(i,i)=sqrt(Q(i,i)-sum(Q(i,1:i-1).^2));
        PP=sum(P(i+1:n-1,1:i-1).*P(i,1:i-1),2);
        QQ=sum(Q(i+1:n-1,1:i-1).*Q(i,1:i-1),2);
        P(i+1:n-1,i)=(P(i+1:n-1,i)-PP(1:n-1-i))/P(i,i);
        Q(i+1:n-1,i)=(Q(i+1:n-1,i)-QQ(1:n-1-i))/Q(i,i);
    end
    P=tril(ones(n-1)).*P;Q=tril(ones(n-1)).*Q;
    Lc=P;
    LL=P.'*Q;
    set(handles.text12,'string','SVD.It takes some time.');drawnow;
    set(handles.text11,'string','(7/12)');drawnow;
    [U,sigma,V]=svd(LL);
    set(handles.text12,'string','Preparing LLL');drawnow;
    set(handles.text11,'string','(8/12)');drawnow;
    LLL=mp(diag(mp(sigma,digitts)),digitts);
    LLL=mp(LLL.^(mp(-1/2,digitts)),digitts);
    set(handles.text12,'string','Preparing T');drawnow;
    set(handles.text11,'string','(9/12)');drawnow;
    T=mp(Lc,digitts)*mp(U,digitts)*mp(diag(LLL),digitts);
    set(handles.text12,'string','Preparing Anew');drawnow;
    set(handles.text11,'string','(10/12)');drawnow;
    Anew=mp(inv(mp(T,digitts)),digitts)*mp(A,digitts)*mp(T,digitts);
    set(handles.text12,'string','Preparing Bnew');drawnow;
    set(handles.text11,'string','(11/12)');drawnow;
    Bnew=mp(inv(mp(T,digitts)),digitts)*mp(B,digitts);
    Cnew=mp(C,digitts)*mp(T,digitts);
    handles.Anew=Anew;
    handles.Bnew=Bnew;
    handles.Cnew=Cnew;
    if handles.popindex==4 || handles.popindex==8
        choos=str2double(get(handles.edit5,'string'));
        handles.mrterms=choos;
        Anew=Anew(1:choos,1:choos);
        Bnew=Bnew(1:choos);
        Cnew=Cnew(1:choos);
        Anew_double=double(Anew);
        [Avec Aeig]=eig(Anew_double);
        Bnew_double=pinv(Avec,1e-14)*Bnew;
        Cnew_double=Cnew*Avec;
        for i=1:choos
            w(i)=Bnew_double(i)*Cnew_double(i);
        end
        p=diag(Aeig);
        approx=zeros(N,1);
        approx=approx+double(X(1));
        if handles.popindex==4 || handles.popindex==5
            v=eval(get(handles.edit8,'string'));
            stepsize=eval(get(handles.edit9,'string'));
            x=(v+stepsize):stepsize:u;
        else
            stepsize=0.002;
            if u>=10
                stepsize=0.2;
            end
            x=stepsize:stepsize:u;
        end
        set(handles.text12,'string','Testing MR');drawnow;
        for i=1:N
            set(handles.text11,'string',[num2str(i/N*100),'%(12/12)']);drawnow;
            for j=1:choos
                approx(i)=approx(i)+w(j)*exp((p(j))*x(i)^2);
            end
        end
        approx=double(approx);
        errors=abs(approx-y);
        handles.errors=max(abs(errors));
        axes(handles.axes1);
        plot(x,errors);
        xlabel('x');ylabel('Error');
        set(handles.text8,'string',num2str(max(abs(approx-y))));drawnow;
        handles.pp=double(p);
        handles.ww=double(w);
        set(handles.text12,'string','Clear!Please save.');drawnow;
        set(handles.text11,'string','100%');drawnow;
        set(handles.pushbutton4,'Visible','On');
        set(handles.text15,'Visible','On');
        set(handles.edit6,'Visible','On');
        guidata(hObject, handles);
    else
        if handles.popindex==5
            v=eval(get(handles.edit8,'string'));
            stepsize=eval(get(handles.edit9,'string'));
            x=(v+stepsize):stepsize:u;
        else
            stepsize=0.002;
            if u>=10
                stepsize=0.2;
            end
            x=stepsize:stepsize:u;
        end
        errs=zeros(1,n);
        xx=x.^2;
        for choos=1:n-1
            set(handles.text12,'string','Testing MR terms');drawnow;
            set(handles.text11,'string',[num2str(choos/n*100) '%(11/11)']);drawnow;
            Anews=Anew(1:choos,1:choos);
            Bnews=Bnew(1:choos);
            Cnews=Cnew(1:choos);
            Anew_double=double(Anews);
            [Avec Aeig]=eig(Anew_double);
            Bnew_double=pinv(Avec,1e-14)*Bnews;
            Cnew_double=Cnews*Avec;
            w=zeros(choos,1);
            for i=1:choos
                w(i)=Bnew_double(i)*Cnew_double(i);
            end
            p=diag(Aeig);
%             approx=zeros(N,1);
%             approx=approx+double(X(1));
%             for i=1:N
%                 for j=1:choos
%                     approx(i)=approx(i)+w(j)*exp((p(j))*x(i)^2);
%                 end
%             end
            approx=sum(w.*exp(p.*xx),1)+double(X(1));
%             errs(choos)=max(abs(approx-y));
            errs(choos)=max(abs(approx-y.'));
        end
        handles.errs=errs;
        axes(handles.axes1);
        semilogy(1:n,errs);
        xlabel('The number of terms');ylabel('Maximum error');
        set(handles.text12,'string','Clear!');drawnow;
        set(handles.text11,'string','100%');drawnow;
        set(handles.text13,'Visible','On');
        set(handles.edit5,'Visible','On');
        set(handles.pushbutton4,'Visible','On');
        set(handles.text15,'Visible','On');
        set(handles.edit6,'Visible','On');
        guidata(hObject, handles); 
    end
elseif handles.popindex==3 || handles.popindex==7
    acc=handles.accindex;
    x=str2double(get(handles.edit3,'string'));handles.xmax=x;
    n=str2double(get(handles.edit1,'string'));handles.n=n;
    if handles.popindex==3
        v=eval(get(handles.edit8,'string'));handles.v=v;
        stepsize=eval(get(handles.edit9,'string'));handles.stepsize=stepsize;
    else
        v=0;stepsize=0.002;
        handles.xmin=v;handles.stepsize=stepsize;
    end
    s=get(handles.edit4,'string');handles.s=s;handles.records=[handles.records,{s}];
    set(handles.popupmenu7,'String',handles.records);
    scstart=eval(get(handles.edit13,'string'));
    scend=eval(get(handles.edit14,'string'));
    scstep=eval(get(handles.edit15,'string'));
    scale=scstart:scstep:scend;
    handles.err=mexptexting2(x,n,scale,s,v,stepsize,acc);
    axes(handles.axes1);
    semilogy(scale,handles.err);
    xlabel('scale/n');ylabel('Error');
    err=min(handles.err);
    whereerr=find(handles.err==err);
    set(handles.text8,'string',num2str(err));drawnow;
    set(handles.text11,'string',['scale/n=' num2str(scale(whereerr(1)))]);drawnow;
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function text8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.popindex=get(handles.popupmenu1,'Value');
if handles.popindex==2
    set(handles.pushbutton1,'Visible','On');
    set(handles.pushbutton3,'Visible','On');
    set(handles.pushbutton4,'Visible','Off');
    set(handles.text3,'Visible','On');
    set(handles.text4,'Visible','On');
    set(handles.text5,'Visible','On');
    set(handles.edit1,'Visible','On');
    set(handles.edit2,'Visible','On');
    set(handles.edit3,'Visible','On');
    set(handles.edit5,'Visible','Off');
    set(handles.text7,'Visible','On');
    set(handles.text8,'Visible','On');
    set(handles.text11,'Visible','Off');
    set(handles.text12,'Visible','Off');
    set(handles.text13,'Visible','Off');
    set(handles.text15,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.edit7,'Visible','Off');
    set(handles.edit8,'Visible','On');
    set(handles.edit9,'Visible','On');
    set(handles.text19,'Visible','On');
    set(handles.text20,'Visible','On');
    set(handles.text23,'Visible','Off');
    set(handles.text24,'Visible','Off');
    set(handles.text25,'Visible','Off');
    set(handles.radiobutton2,'Visible','Off');
    set(handles.text26,'Visible','On');
    set(handles.popupmenu4,'Visible','On');
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
    set(handles.pushbutton6,'Visible','On');
    set(handles.edit13,'Visible','Off');
    set(handles.edit14,'Visible','Off');
    set(handles.edit15,'Visible','Off');
    set(handles.text28,'Visible','Off');
    set(handles.text29,'Visible','Off');
    set(handles.text30,'Visible','Off');
elseif handles.popindex==5
    set(handles.pushbutton1,'Visible','On');
    set(handles.pushbutton3,'Visible','On');
    set(handles.pushbutton4,'Visible','Off');
    set(handles.text3,'Visible','On');
    set(handles.text4,'Visible','On');
    set(handles.text5,'Visible','On');
    set(handles.edit1,'Visible','On');
    set(handles.edit2,'Visible','On');
    set(handles.edit3,'Visible','On');
    set(handles.edit5,'Visible','Off');
    set(handles.text7,'Visible','On');
    set(handles.text8,'Visible','On');
    set(handles.text11,'Visible','On');
    set(handles.text12,'Visible','On');
    set(handles.text13,'Visible','Off');
    set(handles.text15,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.text16,'Visible','On');
    set(handles.edit7,'Visible','On');
    set(handles.edit8,'Visible','On');
    set(handles.edit9,'Visible','On');
    set(handles.text19,'Visible','On');
    set(handles.text20,'Visible','On');
    set(handles.radiobutton2,'Visible','On');
    set(handles.text23,'Visible','Off');
    set(handles.text24,'Visible','Off');
    set(handles.text25,'Visible','Off');
    set(handles.text26,'Visible','On');
    set(handles.popupmenu4,'Visible','On');
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
    set(handles.pushbutton6,'Visible','On');
    set(handles.edit13,'Visible','Off');
    set(handles.edit14,'Visible','Off');
    set(handles.edit15,'Visible','Off');
    set(handles.text28,'Visible','Off');
    set(handles.text29,'Visible','Off');
    set(handles.text30,'Visible','Off');
elseif handles.popindex==3
    set(handles.pushbutton1,'Visible','On');
    set(handles.pushbutton3,'Visible','On');
    set(handles.pushbutton4,'Visible','Off');
    set(handles.text3,'Visible','On');
    set(handles.text4,'Visible','Off');
    set(handles.text5,'Visible','On');
    set(handles.edit1,'Visible','On');
    set(handles.edit2,'Visible','Off');
    set(handles.edit3,'Visible','On');
    set(handles.edit5,'Visible','Off');
    set(handles.text7,'Visible','On');
    set(handles.text8,'Visible','On');
    set(handles.text11,'Visible','On');
    set(handles.text12,'Visible','On');
    set(handles.text13,'Visible','Off');
    set(handles.text15,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.edit7,'Visible','Off');
    set(handles.edit8,'Visible','On');
    set(handles.edit9,'Visible','On');
    set(handles.text19,'Visible','On');
    set(handles.text20,'Visible','On');
    set(handles.radiobutton2,'Visible','Off');
    set(handles.text23,'Visible','Off');
    set(handles.text24,'Visible','Off');
    set(handles.text25,'Visible','Off');
    set(handles.text26,'Visible','On');
    set(handles.popupmenu4,'Visible','On');
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
    set(handles.pushbutton6,'Visible','On');
    set(handles.edit13,'Visible','On');
    set(handles.edit14,'Visible','On');
    set(handles.edit15,'Visible','On');
    set(handles.text28,'Visible','On');
    set(handles.text29,'Visible','On');
    set(handles.text30,'Visible','On');
elseif handles.popindex==4
    set(handles.pushbutton1,'Visible','On');
    set(handles.pushbutton3,'Visible','On');
    set(handles.pushbutton4,'Visible','Off');
    set(handles.text3,'Visible','On');
    set(handles.text4,'Visible','On');
    set(handles.text5,'Visible','On');
    set(handles.text7,'Visible','On');
    set(handles.text8,'Visible','On');
    set(handles.edit1,'Visible','On');
    set(handles.edit2,'Visible','On');
    set(handles.edit3,'Visible','On');
    set(handles.edit5,'Visible','On');
    set(handles.text11,'Visible','On');
    set(handles.text12,'Visible','On');
    set(handles.text13,'Visible','On');
    set(handles.text15,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.text16,'Visible','On');
    set(handles.edit7,'Visible','On');
    set(handles.edit8,'Visible','On');
    set(handles.edit9,'Visible','On');
    set(handles.text19,'Visible','On');
    set(handles.text20,'Visible','On');
    set(handles.radiobutton2,'Visible','On');
    set(handles.text23,'Visible','Off');
    set(handles.text24,'Visible','Off');
    set(handles.text25,'Visible','Off');
    set(handles.text26,'Visible','On');
    set(handles.popupmenu4,'Visible','On');
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
    set(handles.pushbutton6,'Visible','On');
    set(handles.edit13,'Visible','Off');
    set(handles.edit14,'Visible','Off');
    set(handles.edit15,'Visible','Off');
    set(handles.text28,'Visible','Off');
    set(handles.text29,'Visible','Off');
    set(handles.text30,'Visible','Off');
elseif handles.popindex==1  
    set(handles.pushbutton1,'Visible','Off');
    set(handles.pushbutton3,'Visible','Off');
    set(handles.pushbutton4,'Visible','Off');
    set(handles.text3,'Visible','Off');
    set(handles.text4,'Visible','Off');
    set(handles.text5,'Visible','Off');
    set(handles.edit1,'Visible','Off');
    set(handles.edit2,'Visible','Off');
    set(handles.edit3,'Visible','Off');
    set(handles.edit5,'Visible','Off');
    set(handles.text7,'Visible','Off');
    set(handles.text8,'Visible','Off');
    set(handles.text11,'Visible','Off');
    set(handles.text12,'Visible','Off');
    set(handles.text13,'Visible','Off');
    set(handles.text15,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.edit7,'Visible','Off');
    set(handles.edit8,'Visible','Off');
    set(handles.edit9,'Visible','Off');
    set(handles.text19,'Visible','Off');
    set(handles.text20,'Visible','Off');
    set(handles.radiobutton2,'Visible','Off');
    set(handles.text23,'Visible','Off');
    set(handles.text24,'Visible','Off');
    set(handles.text25,'Visible','Off');
    set(handles.text26,'Visible','Off');
    set(handles.popupmenu4,'Visible','Off');
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
    set(handles.pushbutton6,'Visible','On');
    set(handles.edit13,'Visible','Off');
    set(handles.edit14,'Visible','Off');
    set(handles.edit15,'Visible','Off');
    set(handles.text28,'Visible','Off');
    set(handles.text29,'Visible','Off');
    set(handles.text30,'Visible','Off');
elseif handles.popindex==6
    set(handles.pushbutton1,'Visible','On');
    set(handles.pushbutton3,'Visible','On');
    set(handles.pushbutton4,'Visible','Off');
    set(handles.text3,'Visible','On');
    set(handles.text4,'Visible','On');
    set(handles.text5,'Visible','On');
    set(handles.edit1,'Visible','On');
    set(handles.edit2,'Visible','On');
    set(handles.edit3,'Visible','On');
    set(handles.edit5,'Visible','Off');
    set(handles.text7,'Visible','On');
    set(handles.text8,'Visible','On');
    set(handles.text11,'Visible','Off');
    set(handles.text12,'Visible','Off');
    set(handles.text13,'Visible','Off');
    set(handles.text15,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.edit7,'Visible','Off');
    set(handles.edit8,'Visible','Off');
    set(handles.edit9,'Visible','Off');
    set(handles.text19,'Visible','On');
    set(handles.text20,'Visible','On');
    set(handles.text23,'Visible','Off');
    set(handles.text24,'Visible','On');
    set(handles.text25,'Visible','On');
    set(handles.radiobutton2,'Visible','Off');
    set(handles.text26,'Visible','On');
    set(handles.popupmenu4,'Visible','On');
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
    set(handles.pushbutton6,'Visible','On');
    set(handles.edit13,'Visible','Off');
    set(handles.edit14,'Visible','Off');
    set(handles.edit15,'Visible','Off');
    set(handles.text28,'Visible','Off');
    set(handles.text29,'Visible','Off');
    set(handles.text30,'Visible','Off');
elseif handles.popindex==8
    set(handles.pushbutton1,'Visible','On');
    set(handles.pushbutton3,'Visible','On');
    set(handles.pushbutton4,'Visible','Off');
    set(handles.text3,'Visible','On');
    set(handles.text4,'Visible','On');
    set(handles.text5,'Visible','On');
    set(handles.text7,'Visible','On');
    set(handles.text8,'Visible','On');
    set(handles.edit1,'Visible','On');
    set(handles.edit2,'Visible','On');
    set(handles.edit3,'Visible','On');
    set(handles.edit5,'Visible','On');
    set(handles.text11,'Visible','On');
    set(handles.text12,'Visible','On');
    set(handles.text13,'Visible','On');
    set(handles.text15,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.text16,'Visible','On');
    set(handles.edit7,'Visible','Off');
    set(handles.edit8,'Visible','Off');
    set(handles.edit9,'Visible','Off');
    set(handles.text19,'Visible','On');
    set(handles.text20,'Visible','On');
    set(handles.radiobutton2,'Visible','On');
    set(handles.text23,'Visible','On');
    set(handles.text24,'Visible','On');
    set(handles.text25,'Visible','On');
    set(handles.text26,'Visible','On');
    set(handles.popupmenu4,'Visible','On');
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
    set(handles.pushbutton6,'Visible','On');
    set(handles.edit13,'Visible','Off');
    set(handles.edit14,'Visible','Off');
    set(handles.edit15,'Visible','Off');
    set(handles.text28,'Visible','Off');
    set(handles.text29,'Visible','Off');
    set(handles.text30,'Visible','Off');
elseif handles.popindex==9
    set(handles.pushbutton1,'Visible','On');
    set(handles.pushbutton3,'Visible','On');
    set(handles.pushbutton4,'Visible','Off');
    set(handles.text3,'Visible','On');
    set(handles.text4,'Visible','On');
    set(handles.text5,'Visible','On');
    set(handles.edit1,'Visible','On');
    set(handles.edit2,'Visible','On');
    set(handles.edit3,'Visible','On');
    set(handles.edit5,'Visible','Off');
    set(handles.text7,'Visible','On');
    set(handles.text8,'Visible','On');
    set(handles.text11,'Visible','On');
    set(handles.text12,'Visible','On');
    set(handles.text13,'Visible','Off');
    set(handles.text15,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.text16,'Visible','On');
    set(handles.edit7,'Visible','Off');
    set(handles.edit8,'Visible','Off');
    set(handles.edit9,'Visible','Off');
    set(handles.text19,'Visible','On');
    set(handles.text20,'Visible','On');
    set(handles.radiobutton2,'Visible','On');
    set(handles.text23,'Visible','On');
    set(handles.text24,'Visible','On');
    set(handles.text25,'Visible','On');
    set(handles.text26,'Visible','On');
    set(handles.popupmenu4,'Visible','On');
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
    set(handles.pushbutton6,'Visible','On');
    set(handles.edit13,'Visible','Off');
    set(handles.edit14,'Visible','Off');
    set(handles.edit15,'Visible','Off');
    set(handles.text28,'Visible','Off');
    set(handles.text29,'Visible','Off');
    set(handles.text30,'Visible','Off');
elseif handles.popindex==7
    set(handles.pushbutton1,'Visible','On');
    set(handles.pushbutton3,'Visible','On');
    set(handles.pushbutton4,'Visible','Off');
    set(handles.text3,'Visible','On');
    set(handles.text4,'Visible','Off');
    set(handles.text5,'Visible','On');
    set(handles.edit1,'Visible','On');
    set(handles.edit2,'Visible','Off');
    set(handles.edit3,'Visible','On');
    set(handles.edit5,'Visible','Off');
    set(handles.text7,'Visible','On');
    set(handles.text8,'Visible','On');
    set(handles.text11,'Visible','On');
    set(handles.text12,'Visible','On');
    set(handles.text13,'Visible','Off');
    set(handles.text15,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.edit7,'Visible','Off');
    set(handles.edit8,'Visible','Off');
    set(handles.edit9,'Visible','Off');
    set(handles.text19,'Visible','On');
    set(handles.text20,'Visible','On');
    set(handles.radiobutton2,'Visible','Off');
    set(handles.text23,'Visible','Off');
    set(handles.text24,'Visible','On');
    set(handles.text25,'Visible','On');
    set(handles.text26,'Visible','On');
    set(handles.popupmenu4,'Visible','On');
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
    set(handles.pushbutton6,'Visible','On');
    set(handles.edit13,'Visible','On');
    set(handles.edit14,'Visible','On');
    set(handles.edit15,'Visible','On');
    set(handles.text28,'Visible','On');
    set(handles.text29,'Visible','On');
    set(handles.text30,'Visible','On');
end
guidata(hObject, handles); 
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.popindex=get(handles.popupmenu1,'Value');
if handles.popindex==2 || handles.popindex==6
    acc=handles.accindex;
    x=str2double(get(handles.edit3,'string'));
    n=str2double(get(handles.edit1,'string'));
    if handles.popindex==2
        v=eval(get(handles.edit8,'string'));
        stepsize=eval(get(handles.edit9,'string'));
    else
        v=0;stepsize=0.002;
    end
    scale=str2double(get(handles.edit2,'string'));
    s=get(handles.edit4,'string');handles.records=[handles.records,{s}];handles.s=s;
    set(handles.popupmenu7,'String',handles.records);
    set(handles.text8,'string','Now Loading');
    if handles.radioindex==0
        handles.err=exptexting3(x,n,scale,s,v,stepsize,acc);
    else
        s1=s;s2=get(handles.edit10,'string');
        sp=eval(get(handles.edit11,'string'));
        handles.err=exptexting6(x,n,scale,s1,s2,sp,v,stepsize,acc);
    end
    digitts=max(handles.err(3,:));
    set(handles.text8,'string',num2str(digitts));
    axes(handles.axes1);
    plot(handles.err(1,:),handles.err(3,:));
    xlabel('x');ylabel('Error');
    guidata(hObject, handles); 
elseif handles.popindex==4 || handles.popindex==5 || handles.popindex==8 || handles.popindex==9
    u=eval(get(handles.edit3,'string'));handles.xmax=u;
    n=eval(get(handles.edit1,'string'));handles.n=n;
    zz=eval(get(handles.edit2,'string'));handles.scale=zz;
    s=get(handles.edit4,'string');handles.s=[s,'--SOE'];handles.records=[handles.records,{s}];
    set(handles.popupmenu7,'String',handles.records);
    g=@(x) eval(s);
    if handles.popindex==4 || handles.popindex==5
        digitts=eval(get(handles.edit7,'string'));
        handles.digitts=digitts;
    else
        digitts=300;
    end
    mp.Digits(digitts);
    scale=n*zz;
    set(handles.text12,'string','Computing Integrals');drawnow;
    if handles.radioindex==0
        A=zeros(1,2*n);
        if handles.intindex==2
            for i=0:2*n-1
                t=i;
                set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                q=@(x) (2/pi)*cos(t*x).*g(-scale*log((1+cos(x))/2));
                A(i+1)=quadgk(q,0,pi,'AbsTol',handles.accindex,'MaxIntervalCount',1e9,'RelTol',0);
            end
        elseif handles.intindex==3
            w=handles.accindex;
            ww=ceil(sqrt(w));
            d=2;
            S=2*ww/d;
            N=100;
            [b1,b2]=grule(N);
            for i=0:2*n-1
                t=i;
                set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                q=@(x) (2/pi)*cos(t*x).*g(-scale*log((1+cos(x))/2));
                for tt=1:S
                    A(i+1)=A(i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
                end
                A(i+1)=A(i+1)/ww;
            end
        end
    else
        s1=s;s2=get(handles.edit10,'string');
        sp=eval(get(handles.edit11,'string'));
        g1=@(x) eval(s1);g2=@(x) eval(s2);sps=acos(2*exp(-sp/scale)-1);
        A=zeros(1,2*n);
        if handles.intindex==2
            for i=0:2*n-1
                t=i;
                set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                q1=@(x) (2/pi)*cos(t*x).*g1((-scale*log((1+cos(x))/2)));
                q2=@(x) (2/pi)*cos(t*x).*g2((-scale*log((1+cos(x))/2)));
                A(i+1)=quadgk(q1,0,sps,'AbsTol',handles.accindex,'MaxIntervalCount',1e9,'RelTol',0)+quadgk(q2,sps,pi,'AbsTol',1.0e-12,'MaxIntervalCount',1e9,'RelTol',0);
            end
        elseif handles.intindex==3
            w=handles.accindex;
            ww=ceil(sqrt(w));
            d=2;
            S=2*ww/d;
            N=100;
            [b1,b2]=grule(N);
            for i=0:2*n-1
                t=i;
                set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                q1=@(x) (2/pi)*cos(t*x).*g1((-scale*log((1+cos(x))/2)));
                q2=@(x) (2/pi)*cos(t*x).*g2((-scale*log((1+cos(x))/2)));
                for tt=1:S
                    A(i+1)=A(i+1)+q1(((2*tt-1+b1-ww)/ww)*sps/2+sps/2)*b2'*sps/2;
                    A(i+1)=A(i+1)+q2(((2*tt-1+b1-ww)/ww)*(pi-sps)/2+(sps+pi)/2)*b2'*(pi-sps)/2;
                end
                A(i+1)=A(i+1)/ww;
            end
        end
    end
    A=mp(A,digitts);
    maxC=0.00;
    set(handles.text12,'string','Computing Linear Coefficients(VP)');drawnow;
    for r=0:2*n-1
        set(handles.text11,'string',[num2str(r/2/n*100),'%(2/12)']);drawnow;
        if(r==0)
            sum1=mp('0.00');
            for k=1:n
                sum1=mp(sum1+mp(-1)^mp(k)*mp(A(k+1)),digitts);
            end
            sum2=mp('0.00');
            for k=1:n-1
                sum2=mp(sum2+mp(-1)^mp(n+k)*(1-mp(k)/n)*mp(A(n+k+1)),digitts);
            end
            C(r+1)=mp(A(r+1)/2+sum1+sum2,digitts);
        end
        if(1<=r&&r<=n)
            sum1=mp('0.00');
            sum2=mp('0.00');
            for k=r:n
                sum1=mp(sum1+mp((-1)^(k-r))*mp(2*mp(k)/mp(k+r))* mp(factorial(mp(k+r))) / mp(factorial(mp(2*r))*factorial(mp(k-r)))*mp((mp(2)^mp(2*r-1)))*mp(A(k+1)),digitts );
            end
            for k=1:n-1
                sum2=mp(sum2+mp((-1)^(n+k-r))*mp(1-mp(k)/n)*2*mp((n+k)*(mp(1)/(n+k+r)))*(mp(factorial(mp(n+k+r))))/((mp(factorial(mp(2*r))))*(mp(factorial(mp(n+k-r)))))*mp(2^mp(2*r-1))*mp(A(n+k+1)),digitts);
            end
            C(r+1)=mp((sum1)+(sum2),digitts);
        end
        if(r>=n+1)
            sum1=mp('0.00');
            for k=(r-n:n-1)
                sum1=mp(sum1+mp(mp(-1)^mp(n+k-r))*mp(1-mp(k)/n)*2*mp(mp(n+k)/(mp(n+k+r)))*(mp(factorial(mp(n+k+r))))/(mp(factorial(mp(2*r)))*mp(factorial(mp(n+k-r))))*(mp(2^mp(2*r-1)))*mp(A(n+k+1)),digitts);
            end
            C(r+1)=mp(sum1,digitts);
        end
    end
    if handles.popindex==4 || handles.popindex==5
        v=eval(get(handles.edit8,'string'));
        stepsize=eval(get(handles.edit9,'string'));
        handles.xmin=v;handles.stepsize=stepsize;
        x=(v+stepsize):stepsize:u;
    else
        stepsize=0.002;v=0;
        handles.xmin=v;handles.stepsize=stepsize;
        if u>=10
            stepsize=0.2;
        end
        x=stepsize:stepsize:u;
    end
    if handles.radioindex2==0
        NNN=size(x,2);
%         for i=1:NNN
%             i;
%             a(i)=mp('0.00');
%         end
        set(handles.text12,'string','Testing the VP sum');drawnow;
%         for k=0:2*n-1
%             set(handles.text11,'string',[num2str(k/2/n*100),'%(2/11)']);drawnow;
%             for i=1:NNN
%                 y(i)=mp(exp(mp(-mp(k)*(mp(x(i)))/mp(scale))));
%                 a(i)=mp(a(i)+mp(C(k+1))*mp(y(i)),digitts);
%             end 
%         end
        set(handles.text11,'string','(3/12)');drawnow;
        CC=C.';xx=mp(x);KK=mp(0:(2*n-1))/mp(scale);
        a=sum(CC.*exp(-KK'.*xx),1);
        eta=double(max(mp(abs(mp(a-(g(x)))))));
%         eta=double(max(mp(abs(mp(a-(g(x)))))));
        if handles.radioindex==1
            eta1=[ones(1,round((sp-v)/stepsize)) zeros(1,round(NNN-(sp-v)/stepsize))];
            eta2=[zeros(1,round((sp-v)/stepsize)) ones(1,round(NNN-(sp-v)/stepsize))];
            gg=g1(x).*eta1+g2(x).*eta2;
            eta=double(max(mp(abs(mp(a-gg)))));
        end
        set(handles.text8,'string',[num2str(eta) '(VP)']);drawnow;
    end
    X=C;handles.XX=double(X(1));
    %%Model Reduction
    set(handles.text12,'string','Model Reduction');drawnow;
    set(handles.text11,'string','(4/12)');drawnow;
    n=2*n;mp.Digits(digitts);
    if handles.popindex==4 || handles.popindex==5
        v=eval(get(handles.edit8,'string'));
        stepsize=eval(get(handles.edit9,'string'));
        x=(v+stepsize):stepsize:u;
    else
        stepsize=0.002;
        if u>=10
            stepsize=0.2;
        end
        x=stepsize:stepsize:u;
    end
    N=size(x,2);
    y=g(x);
    if handles.radioindex==1
        y=g1(x).*eta1+g2(x).*eta2;
    end
    y=y.';
    A=mp(-mp(diag(1:n-1))/mp(n*zz/2),digitts);
    A2=mp(-mp((1:n-1))/mp(n*zz/2),digitts);
    B=mp(zeros(n-1,1),digitts);
    C=mp(zeros(1,n-1),digitts);
    for i=1:n-1
        if(X(i+1)>0)
            C(i)=mp(sqrt(mp(X(i+1),digitts)),digitts);
            B(i)=mp(sqrt(mp(X(i+1),digitts)),digitts);
        else
            C(i)=mp(sqrt(mp(-X(i+1),digitts)),digitts);
            B(i)=mp(-sqrt(mp(-X(i+1),digitts)),digitts);
        end
    end
    BBB=B*B.';
    CCC=C.'*C;
    AA=A2+A2.';
    set(handles.text12,'string','Preparing P and Q');drawnow;
    set(handles.text11,'string',' ');drawnow;
    P=-BBB./AA;Q=-CCC./AA;
    set(handles.text12,'string','Cholesky');drawnow;
    mp.Digits(digitts);
    for i=1:n-1
        set(handles.text11,'string',[num2str(i/n*100),'%(6/12)']);drawnow;
        P(i,i)=sqrt(P(i,i)-sum(P(i,1:i-1).^2));
        Q(i,i)=sqrt(Q(i,i)-sum(Q(i,1:i-1).^2));
        PP=sum(P(i+1:n-1,1:i-1).*P(i,1:i-1),2);
        QQ=sum(Q(i+1:n-1,1:i-1).*Q(i,1:i-1),2);
        P(i+1:n-1,i)=(P(i+1:n-1,i)-PP(1:n-1-i))/P(i,i);
        Q(i+1:n-1,i)=(Q(i+1:n-1,i)-QQ(1:n-1-i))/Q(i,i);
    end
    P=tril(ones(n-1)).*P;Q=tril(ones(n-1)).*Q;
    Lc=P;
    LL=P.'*Q;
    set(handles.text12,'string','SVD.It takes some time.');drawnow;
    set(handles.text11,'string','(7/12)');drawnow;
    [U,sigma,V]=svd(LL);
    set(handles.text12,'string','Preparing LLL');drawnow;
    set(handles.text11,'string','(8/12)');drawnow;
    LLL=mp(diag(mp(sigma,digitts)),digitts);
    LLL=mp(LLL.^(mp(-1/2,digitts)),digitts);
    set(handles.text12,'string','Preparing T');drawnow;
    set(handles.text11,'string','(9/12)');drawnow;
    T=mp(Lc,digitts)*mp(U,digitts)*mp(diag(LLL),digitts);
    set(handles.text12,'string','Preparing Anew');drawnow;
    set(handles.text11,'string','(10/12)');drawnow;
    Anew=mp(inv(mp(T,digitts)),digitts)*mp(A,digitts)*mp(T,digitts);
    set(handles.text12,'string','Preparing Bnew');drawnow;
    set(handles.text11,'string','(11/12)');drawnow;
    Bnew=mp(inv(mp(T,digitts)),digitts)*mp(B,digitts);
    Cnew=mp(C,digitts)*mp(T,digitts);
    handles.Anew=Anew;
    handles.Bnew=Bnew;
    handles.Cnew=Cnew;
    if handles.popindex==4 || handles.popindex==8
        choos=str2double(get(handles.edit5,'string'));
        handles.mrterms=choos;
        Anew=Anew(1:choos,1:choos);
        Bnew=Bnew(1:choos);
        Cnew=Cnew(1:choos);
        Anew_double=double(Anew);
        [Avec Aeig]=eig(Anew_double);
        Bnew_double=pinv(Avec,1e-14)*Bnew;
        Cnew_double=Cnew*Avec;
        for i=1:choos
            w(i)=Bnew_double(i)*Cnew_double(i);
        end
        p=diag(Aeig);
        approx=zeros(N,1);
        approx=approx+double(X(1));
        if handles.popindex==4 || handles.popindex==5
            v=eval(get(handles.edit8,'string'));
            stepsize=eval(get(handles.edit9,'string'));
            x=(v+stepsize):stepsize:u;
        else
            stepsize=0.002;
            if u>=10
                stepsize=0.2;
            end
            x=stepsize:stepsize:u;
        end
        set(handles.text12,'string','Testing MR');drawnow;
        for i=1:N
            set(handles.text11,'string',[num2str(i/N*100),'%(12/12)']);drawnow;
            for j=1:choos
                approx(i)=approx(i)+w(j)*exp((p(j))*x(i));
            end
        end
        approx=double(approx);
        errors=abs(approx-y);
        axes(handles.axes1);
        plot(x,errors);
        xlabel('x');ylabel('Error');
        handles.errors=max(abs(errors));
        set(handles.text8,'string',num2str(max(abs(approx-y))));drawnow;
        handles.pp=double(p);
        handles.ww=double(w);
        set(handles.text12,'string','Clear!Please save.');drawnow;
        set(handles.text11,'string','100%');drawnow;
        set(handles.pushbutton4,'Visible','On');
        set(handles.text15,'Visible','On');
        set(handles.edit6,'Visible','On');
        guidata(hObject, handles);
    else
        if handles.popindex==5
            v=eval(get(handles.edit8,'string'));
            stepsize=eval(get(handles.edit9,'string'));
            x=(v+stepsize):stepsize:u;
        else
            stepsize=0.002;
            if u>=10
                stepsize=0.2;
            end
            x=stepsize:stepsize:u;
        end
        errs=zeros(1,n);
        for choos=1:n-1
            set(handles.text12,'string','Testing MR terms');drawnow;
            set(handles.text11,'string',[num2str(choos/n*100) '%(11/11)']);drawnow;
            Anews=Anew(1:choos,1:choos);
            Bnews=Bnew(1:choos);
            Cnews=Cnew(1:choos);
            Anew_double=double(Anews);
            [Avec Aeig]=eig(Anew_double);
            Bnew_double=pinv(Avec,1e-14)*Bnews;
            Cnew_double=Cnews*Avec;
            w=zeros(choos,1);
            for i=1:choos
                w(i)=Bnew_double(i)*Cnew_double(i);
            end
            p=diag(Aeig);
%             approx=zeros(N,1);
%             approx=approx+double(X(1));
%             for i=1:N
%                 for j=1:choos
%                     approx(i)=approx(i)+w(j)*exp((p(j))*x(i));
%                 end
%             end
            approx=sum(w.*exp(p.*x),1)+double(X(1));
            errs(choos)=max(abs(approx-y.'));
        end
        handles.errs=errs;
        axes(handles.axes1);
        semilogy(1:n,errs);
        xlabel('The number of terms');ylabel('Maximum Error');
        set(handles.text12,'string','Clear!');drawnow;
        set(handles.text11,'string','100%');drawnow;
        set(handles.text13,'Visible','On');
        set(handles.edit5,'Visible','On');
        set(handles.pushbutton4,'Visible','On');
        set(handles.text15,'Visible','On');
        set(handles.edit6,'Visible','On');
        guidata(hObject, handles); 
    end
elseif handles.popindex==3 || handles.popindex==7
    acc=handles.accindex;
    x=str2double(get(handles.edit3,'string'));handles.xmax=x;
    n=str2double(get(handles.edit1,'string'));handles.n=n;
    if handles.popindex==7
        v=eval(get(handles.edit8,'string'));handles.v=v;
        stepsize=eval(get(handles.edit9,'string'));handles.stepsize=stepsize;
    else
        v=0;stepsize=0.002;
        handles.xmin=v;handles.stepsize=stepsize;
    end
    guidata(hObject, handles);
    s=get(handles.edit4,'string');handles.s=s;handles.records=[handles.records,{s}];
    set(handles.popupmenu7,'String',handles.records);
    scstart=eval(get(handles.edit13,'string'));
    scend=eval(get(handles.edit14,'string'));
    scstep=eval(get(handles.edit15,'string'));
    scale=scstart:scstep:scend;
    handles.err=mexptexting3(x,n,scale,s,v,stepsize,acc);
    axes(handles.axes1);
    semilogy(scale,handles.err);
    xlabel('scale/n');ylabel('Error');
    err=min(handles.err);
    whereerr=find(handles.err==err);
    set(handles.text8,'string',num2str(err));drawnow;
    set(handles.text11,'string',['scale/n=' num2str(scale(whereerr(1)))]);drawnow;
end
%(besselh(0,1,50*(x+0.05)).*exp(-50j*(x+0.05)))/4

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Xs=handles.XX;
ps=handles.pp;
ws=handles.ww;
if handles.popindex==5 || handles.popindex==9
    choos=str2double(get(handles.edit5,'string'));
    handles.mrterms=choos;
    Anew=handles.Anew;
    Bnew=handles.Bnew;
    Cnew=handles.Cnew;
    Anews=Anew(1:choos,1:choos);
    Bnews=Bnew(1:choos);
    Cnews=Cnew(1:choos);
    Anew_double=double(Anews);
    [Avec Aeig]=eig(Anew_double);
    Bnew_double=pinv(Avec,1e-14)*Bnews;
    Cnew_double=Cnews*Avec;
    for i=1:choos
        ws(i)=Bnew_double(i)*Cnew_double(i);
    end
    ps=diag(Aeig);
    errs=handles.errs;
    handles.errors=errs(choos);
end
pathname=[get(handles.edit6,'string'),'/parameter',datestr(now,30)];
mkdir(pathname);
save([pathname,'/X'],'Xs');
save([pathname,'/p'],'ps');
save([pathname,'/w'],'ws');
fd=fopen([pathname,'/readme.txt'],'w');
t0='function=';
t1='n=';t2='scale/n=';t3='digits=';
t4='step=';t5='xmin=';t6='xmax=';
t7='MRterms=';t8='Maxerror=';
fprintf(fd,'%s',t0);fprintf(fd,'%s',handles.s);fprintf(fd,'\r\n');
fprintf(fd,'%s',t1);fprintf(fd,'%g',handles.n);fprintf(fd,'\r\n');
fprintf(fd,'%s',t2);fprintf(fd,'%g',handles.scale);fprintf(fd,'\r\n');
fprintf(fd,'%s',t3);fprintf(fd,'%g',handles.digitss);fprintf(fd,'\r\n');
fprintf(fd,'%s',t4);fprintf(fd,'%g',handles.stepsize);fprintf(fd,'\r\n');
fprintf(fd,'%s',t5);fprintf(fd,'%g',handles.xmin);fprintf(fd,'\r\n');
fprintf(fd,'%s',t6);fprintf(fd,'%g',handles.xmax);fprintf(fd,'\r\n');
fprintf(fd,'%s',t7);fprintf(fd,'%g',handles.mrterms);fprintf(fd,'\r\n');
fprintf(fd,'%s',t8);fprintf(fd,'%g',handles.errors);fprintf(fd,'\r\n');
fclose(fd);
guidata(hObject, handles); 

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.radioindex=get(handles.radiobutton1,'Value');
if handles.radioindex==1
    set(handles.text22,'Visible','On');
    set(handles.edit10,'Visible','On');
    set(handles.edit11,'Visible','On');
else
    set(handles.text22,'Visible','Off');
    set(handles.edit10,'Visible','Off');
    set(handles.edit11,'Visible','Off');
end
guidata(hObject, handles); 
% Hint: get(hObject,'Value') returns toggle state of radiobutton1



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.radioindex2=get(handles.radiobutton2,'Value');
guidata(hObject, handles); 
% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes during object creation, after setting all properties.
function text23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function ys=exptexting2(u,n,scale,s,v,stepsize,acc)
%SOG
g=@(x) eval(s);
scale=n*scale;
A=zeros(1,2*n);
if acc<1
    for i=0:2*n-1
        t=i;
        q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
        A(i+1)=quadgk(q,0,pi,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0);
    end
else
    w=acc;
    ww=ceil(sqrt(w));
    d=2;
    S=2*ww/d;
    N=100;
    [b1,b2]=grule(N);
    for i=0:2*n-1
        t=i;
        q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
        for tt=1:S
            A(i+1)=A(i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
        end
        A(i+1)=A(i+1)/ww;
    end
end
A(1)=A(1)/2;
u=(v+stepsize):stepsize:u;
y=2*exp(-(u.^2)/scale)-1;
for i=1:size(u,2)
    T=Chebyshev(y(i),2*n-1);
    K(i)=0;
    for k=0:n
        K(i)=K(i)+A(k+1)*T(k+1);
    end
    for k=1:n-1
        K(i)=K(i)+(1-k/n)*A(k+n+1)*T(k+n+1);
    end
end
z=g(u);
y=abs(z-K);
ys=[u;K;y];

function ys=exptexting3(u,n,scale,s,v,stepsize,acc)
g=@(x) eval(s);
scale=n*scale;
A=zeros(1,2*n);
if acc<1
    for i=0:2*n-1
        t=i;
        q=@(x) (2/pi)*cos(t*x).*g((-scale*log((1+cos(x))/2)));
        A(i+1)=quadgk(q,0,pi,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0);
    end
else
    w=acc;
    ww=ceil(sqrt(w));
    d=2;
    S=2*ww/d;
    N=100;
    [b1,b2]=grule(N);
    for i=0:2*n-1
        t=i;
        q=@(x) (2/pi)*cos(t*x).*g((-scale*log((1+cos(x))/2)));
        for tt=1:S
            A(i+1)=A(i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
        end
        A(i+1)=A(i+1)/ww;
    end
end
A(1)=A(1)/2;
u=(v+stepsize):stepsize:u;
y=2*exp(-(u)/scale)-1;
for i=1:size(u,2)
    T=Chebyshev(y(i),2*n-1);
    K(i)=0;
    for k=0:n
        K(i)=K(i)+A(k+1)*T(k+1);
    end
    for k=1:n-1
        K(i)=K(i)+(1-k/n)*A(k+n+1)*T(k+n+1);
    end
end
z=g(u);
y=abs(z-K);
ys=[u;K;y];

function ys=exptexting5(u,n,scale,s1,s2,sp,v,stepsize,acc)
g1=@(x) eval(s1);
g2=@(x) eval(s2);
scale=n*scale;sps=acos(2*exp(-sp*sp/scale)-1);
A=zeros(1,2*n);
if acc<1
    for i=0:2*n-1
        t=i;
        q1=@(x) (2/pi)*cos(t*x).*g1(sqrt(-scale*log((1+cos(x))/2)));
        q2=@(x) (2/pi)*cos(t*x).*g2(sqrt(-scale*log((1+cos(x))/2)));
        A(i+1)=quadgk(q1,0,sps,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0)+quadgk(q2,sps,pi,'AbsTol',1.0e-12,'MaxIntervalCount',1e9,'RelTol',0);
    end
else
    w=acc;
    ww=ceil(sqrt(w));
    d=2;
    S=2*ww/d;
    N=100;
    [b1,b2]=grule(N);
    for i=0:2*n-1
        t=i;
        set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
        q1=@(x) (2/pi)*cos(t*x).*g1(sqrt(-scale*log((1+cos(x))/2)));
        q2=@(x) (2/pi)*cos(t*x).*g2(sqrt(-scale*log((1+cos(x))/2)));
        for tt=1:S
            A(i+1)=A(i+1)+q1(((2*tt-1+b1-ww)/ww)*sps/2+sps/2)*b2'*sps/2;
            A(i+1)=A(i+1)+q2(((2*tt-1+b1-ww)/ww)*(pi-sps)/2+(sps+pi)/2)*b2'*(pi-sps)/2;
        end
        A(i+1)=A(i+1)/ww;
    end
end
A(1)=A(1)/2;
u=(v+stepsize):stepsize:u;
NNN=size(u,2);
y=2*exp(-(u.^2)/scale)-1;
for i=1:size(u,2)
    T=Chebyshev(y(i),2*n-1);
    K(i)=0;
    for k=0:n
        K(i)=K(i)+A(k+1)*T(k+1);
    end
    for k=1:n-1
        K(i)=K(i)+(1-k/n)*A(k+n+1)*T(k+n+1);
    end
end
eta1=[ones(1,(sp-v)/stepsize) zeros(1,NNN-(sp-v)/stepsize)];
eta2=[zeros(1,(sp-v)/stepsize) ones(1,NNN-(sp-v)/stepsize)];
z=g1(u).*eta1+g2(u).*eta2;
y=abs(z-K);
ys=[u;K;y];

function ys=exptexting6(u,n,scale,s1,s2,sp,v,stepsize,acc)
g1=@(x) eval(s1);
g2=@(x) eval(s2);
scale=n*scale;sps=acos(2*exp(-sp/scale)-1);
A=zeros(1,2*n);
if acc<1
    for i=0:2*n-1
        t=i;
        q1=@(x) (2/pi)*cos(t*x).*g1((-scale*log((1+cos(x))/2)));
        q2=@(x) (2/pi)*cos(t*x).*g2((-scale*log((1+cos(x))/2)));
        A(i+1)=quadgk(q1,0,sps,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0)+quadgk(q2,sps,pi,'AbsTol',1.0e-12,'MaxIntervalCount',1e9,'RelTol',0);
    end
else
    w=acc;
    ww=ceil(sqrt(w));
    d=2;
    S=2*ww/d;
    N=100;
    [b1,b2]=grule(N);
    for i=0:2*n-1
        t=i;
        set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
        q1=@(x) (2/pi)*cos(t*x).*g1((-scale*log((1+cos(x))/2)));
        q2=@(x) (2/pi)*cos(t*x).*g2((-scale*log((1+cos(x))/2)));
        for tt=1:S
            A(i+1)=A(i+1)+q1(((2*tt-1+b1-ww)/ww)*sps/2+sps/2)*b2'*sps/2;
            A(i+1)=A(i+1)+q2(((2*tt-1+b1-ww)/ww)*(pi-sps)/2+(sps+pi)/2)*b2'*(pi-sps)/2;
        end
        A(i+1)=A(i+1)/ww;
    end
end
A(1)=A(1)/2;
u=(v+stepsize):stepsize:u;
NNN=size(u,2);
y=2*exp(-(u)/scale)-1;
for i=1:size(u,2)
    T=Chebyshev(y(i),2*n-1);
    K(i)=0;
    for k=0:n
        K(i)=K(i)+A(k+1)*T(k+1);
    end
    for k=1:n-1
        K(i)=K(i)+(1-k/n)*A(k+n+1)*T(k+n+1);
    end
end
eta1=[ones(1,(sp-v)/stepsize) zeros(1,NNN-(sp-v)/stepsize)];
eta2=[zeros(1,(sp-v)/stepsize) ones(1,NNN-(sp-v)/stepsize)];
z=g1(u).*eta1+g2(u).*eta2;
y=abs(z-K);
ys=[u;K;y];

function a = Chebyshev(x,n)
%代入n，计算Tn
xsize=size(x,2);
a=zeros(xsize,n+1);
a(:,1)=1;
a(:,2)=x.';
for i=3:n+1
    a(:,i)=2*x.'.*a(:,i-1)-a(:,i-2);
end



% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.intindex=get(handles.popupmenu4,'Value');
p=get(handles.popupmenu4,'Value');
if p==2
    set(handles.popupmenu5,'Visible','On');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','On');
elseif p==3
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','On');
    set(handles.text27,'Visible','On');
elseif p==1
    set(handles.popupmenu5,'Visible','Off');
    set(handles.popupmenu6,'Visible','Off');
    set(handles.text27,'Visible','Off');
end
guidata(hObject, handles); 
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [bp,wf]=grule(n)
% [bp,wf]=grule(n)
%  This function computes Gauss base points and weight factors
%  using the algorithm given by Davis and Rabinowitz in 'Methods
%  of Numerical Integration', page 365, Academic Press, 1975.
bp=zeros(n,1); wf=bp; iter=2; m=fix((n+1)/2); e1=n*(n+1);
mm=4*m-1; t=(pi/(4*n+2))*(3:4:mm); nn=(1-(1-1/n)/(8*n*n));
xo=nn*cos(t);
for j=1:iter
   pkm1=1; pk=xo;
   for k=2:n
      t1=xo.*pk; pkp1=t1-pkm1-(t1-pkm1)/k+t1;
      pkm1=pk; pk=pkp1;
   end
   den=1.-xo.*xo; d1=n*(pkm1-xo.*pk); dpn=d1./den;
   d2pn=(2.*xo.*dpn-e1.*pk)./den;
   d3pn=(4*xo.*d2pn+(2-e1).*dpn)./den;
   d4pn=(6*xo.*d3pn+(6-e1).*d2pn)./den;
   u=pk./dpn; v=d2pn./dpn;
   h=-u.*(1+(.5*u).*(v+u.*(v.*v-u.*d3pn./(3*dpn))));
   p=pk+h.*(dpn+(.5*h).*(d2pn+(h/3).*(d3pn+.25*h.*d4pn)));
   dp=dpn+h.*(d2pn+(.5*h).*(d3pn+h.*d4pn/3));
   h=h-p./dp; xo=xo+h;
end
bp=-xo-h;
fx=d1-h.*e1.*(pk+(h/2).*(dpn+(h/3).*(...
    d2pn+(h/4).*(d3pn+(.2*h).*d4pn))));
wf=2*(1-bp.^2)./(fx.*fx);
if (m+m) > n, bp(m)=0; end
if ~((m+m) == n), m=m-1; end
jj=1:m; n1j=(n+1-jj); bp(n1j)=-bp(jj); wf(n1j)=wf(jj);


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p=get(handles.popupmenu5,'Value');
if p==2
    handles.accindex=1e-8;
elseif p==3
    handles.accindex=1e-10;
elseif p==4
    handles.accindex=1e-12;
elseif p==5
    handles.accindex=1e-14;
end
guidata(hObject, handles); 
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p=get(handles.popupmenu6,'Value');
if p==2
    handles.accindex=5e6;
elseif p==3
    handles.accindex=1e7;
elseif p==4
    handles.accindex=5e7;
elseif p==5
    handles.accindex=1e8;
end
guidata(hObject, handles); 
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.popindex=1;
handles.err=[];
handles.XX=0;
handles.pp=[];
handles.ww=[];
handles.Anew=[];
handles.Bnew=[];
handles.Cnew=[];
handles.radioindex=0;
handles.radioindex2=0;
handles.digitss=300;
handles.xmax=0;
handles.xmin=0;
handles.n=0;
handles.stepsize=0;
handles.sclae=0;
handles.errors=0;
handles.errs=[];
handles.s=[];
handles.accindex=1;
handles.mrterms=0;
handles.intindex=2;
set(handles.pushbutton1,'Visible','Off');
set(handles.pushbutton3,'Visible','Off');
set(handles.pushbutton4,'Visible','Off');
set(handles.pushbutton6,'Visible','Off');
set(handles.text3,'Visible','Off');
set(handles.text4,'Visible','Off');
set(handles.text5,'Visible','Off');
set(handles.edit1,'Visible','Off');
set(handles.edit2,'Visible','Off');
set(handles.edit3,'Visible','Off');
set(handles.text7,'Visible','Off');
set(handles.text8,'Visible','Off');
set(handles.text11,'Visible','Off');
set(handles.text13,'Visible','Off');
set(handles.edit5,'Visible','Off');
set(handles.text15,'Visible','Off');
set(handles.edit6,'Visible','Off');
set(handles.text16,'Visible','Off');
set(handles.text19,'Visible','Off');
set(handles.text20,'Visible','Off');
set(handles.edit7,'Visible','Off');
set(handles.edit8,'Visible','Off');
set(handles.edit9,'Visible','Off');
set(handles.edit10,'Visible','Off');
set(handles.edit11,'Visible','Off');
set(handles.text22,'Visible','Off');
set(handles.text23,'Visible','Off');
set(handles.text24,'Visible','Off');
set(handles.text25,'Visible','Off');
set(handles.text26,'Visible','Off');
set(handles.popupmenu4,'Visible','Off');
set(handles.popupmenu5,'Visible','Off');
set(handles.popupmenu6,'Visible','Off');
set(handles.text27,'Visible','Off');
set(handles.radiobutton2,'Visible','Off');
% Update handles structure
guidata(hObject, handles);



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function err=mexptexting2(u,n,scales,s,v,stepsize,acc)
%SOG
g=@(x) eval(s);
scales=n*scales;
sc=size(scales,2);
A=zeros(sc,2*n);
if acc<1
    for j=1:sc
        scale=scales(j);
        for i=0:2*n-1
            t=i;
            q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
            A(j,i+1)=quadgk(q,0,pi,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0);
        end
    end
else
     w=acc;
    ww=ceil(sqrt(w));
    d=2;
    S=2*ww/d;
    N=100;
    [b1,b2]=grule(N);
    for j=1:sc
        scale=scales(j);
        for i=0:2*n-1
            t=i;
            q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
            for tt=1:S
                A(j,i+1)=A(j,i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
            end
            A(j,i+1)=A(j,i+1)/ww;
        end
    end
end
A(:,1)=A(:,1)/2;
u=(v+stepsize):stepsize:u;
scales=scales.';
y=2*exp(-(u.^2)./scales)-1;
K=zeros(sc,size(u,2));
for i=1:size(u,2)
    T=Chebyshev(y(:,i).',2*n-1);
    for k=0:n
        K(:,i)=K(:,i)+A(:,k+1).*T(:,k+1);
    end
    for k=1:n-1
        K(:,i)=K(:,i)+(1-k/n)*A(:,k+n+1).*T(:,k+n+1);
    end
end
z=g(u).*ones(sc,1);
y=abs(z-K);
err=max(y.');

function err=mexptexting3(u,n,scales,s,v,stepsize,acc)
%SOG
g=@(x) eval(s);
scales=n*scales;
sc=size(scales,2);
A=zeros(sc,2*n);
if acc<1
    for j=1:sc
        scale=scales(j);
        for i=0:2*n-1
            t=i;
            q=@(x) (2/pi)*cos(t*x).*g((-scale*log((1+cos(x))/2)));
            A(j,i+1)=quadgk(q,0,pi,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0);
        end
    end
else
     w=acc;
    ww=ceil(sqrt(w));
    d=2;
    S=2*ww/d;
    N=100;
    [b1,b2]=grule(N);
    for j=1:sc
        scale=scales(j);
        for i=0:2*n-1
            t=i;
            q=@(x) (2/pi)*cos(t*x).*g((-scale*log((1+cos(x))/2)));
            for tt=1:S
                A(j,i+1)=A(j,i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
            end
            A(j,i+1)=A(j,i+1)/ww;
        end
    end
end
A(:,1)=A(:,1)/2;
u=(v+stepsize):stepsize:u;
scales=scales.';
y=2*exp(-(u)./scales)-1;
K=zeros(sc,size(u,2));
for i=1:size(u,2)
    T=Chebyshev(y(:,i).',2*n-1);
    for k=0:n
        K(:,i)=K(:,i)+A(:,k+1).*T(:,k+1);
    end
    for k=1:n-1
        K(:,i)=K(:,i)+(1-k/n)*A(:,k+n+1).*T(:,k+n+1);
    end
end
z=g(u).*ones(sc,1);
y=abs(z-K);
err=max(y.');


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index2=get(handles.popupmenu7,'Value');
handles.s=handles.records{1,index2};
set(handles.edit4,'String',handles.s);
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
