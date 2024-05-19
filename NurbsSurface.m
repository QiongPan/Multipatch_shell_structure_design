function multiPatSurf = NurbsSurface(type)

%%%%%% stiffness matrix for each patch, computed elementwise on Guass integration points %%%%%%
% input:
%   type            - the label of surface
%                   - 1-2-Hyperbolic, 4-torus, 5-table holder, 6-pipes, 3-bird's nest
% output:
%   multiPatSurf    - the geometric information of the mid-surface
%                   - a cell with dim = num_pat * 1
%                   - multiPatSurf(k,:) = {p,q,m,n,U,V,Ctrlpts}

switch type 
    case 1  % Hyperbolic 1
        Npat = 1;
        multiPatSurf = cell(Npat,7); % p,q,m,n,U,V,Ctrlpts.
        
        %% patch 1
        p = 2;  q = 2;
        m = 2;  n = 2;
        % contorl point
        Points = [
            -12, -12, 0
            -12, 0, 8
            -12, 12, 0
            0, -12, 8
            0, 0, 16
            0, 12, 8
            12, -12, 0
            12, 0, 8
            12, 12, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            1
            1
            1
            1
            ];
        % knots
        U = [
            0
            0
            1
            1
            ];
        
        V = [
            0
            0
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(1,:) = {p,q,m,n,U,V,Ctrlpts} ;
    case 2  % Hyperbolic 4
        Npat = 4;
        multiPatSurf = cell(Npat,7); % p,q,m,n,U,V,Ctrlpts.
        
        %% patch 1
        p = 2;  q = 2;
        m = 2;  n = 2;
        % contorl point
        Points = [
            -12, -12, 0
            -12, -4.8, 4.8
            -12, 2.4, 3.84
            -7.8, -12, 2.8
            -7.8, -4.8, 7.6
            -7.8, 2.4, 6.64
            -3.6, -12, 3.64
            -3.6, -4.8, 8.44
            -3.6, 2.4, 7.48   ];
        % weights
        W = [
            1
            1
            1
            1
            1
            1
            1
            1
            1
            ];
        % knots
        U = [
            0
            0
            1
            1
            ];
        
        V = [
            0
            0
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(1,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 2
        p = 2;  q = 2;
        m = 2;  n = 2;
        % contorl point
        Points = [
            -3.6, -12, 3.64
            -3.6, -4.8 8.44
            -3.6, 2.4, 7.48
            4.2, -12, 5.2
            4.2, -4.8, 10
            4.2, 2.4, 9.04
            12, -12, 0
            12, -4.8 4.8
            12, 2.4, 3.84  ];
        % weights
        W = [
            1
            1
            1
            1
            1
            1
            1
            1
            1
            ];
        % knots
        U = [
            0
            0
            1
            1
            ];
        
        V = [
            0
            0
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(2,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 3
        p = 2;  q = 2;
        m = 2;  n = 2;
        % contorl point
        Points = [
            -3.6, 2.4, 7.48
            -3.6, 7.2, 6.84
            -3.6, 12, 3.64
            4.2, 2.4, 9.04
            4.2, 7.2, 8.4
            4.2, 12, 5.2
            12, 2.4, 3.84
            12, 7.2, 3.2
            12, 12, 0  ];
        % weights
        W = [
            1
            1
            1
            1
            1
            1
            1
            1
            1
            ];
        % knots
        U = [
            0
            0
            1
            1
            ];
        
        V = [
            0
            0
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(3,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 4
        p = 2;  q = 2;
        m = 2;  n = 2;
        % contorl point
        Points = [
            -12, 2.4, 3.84
            -12, 7.2, 3.2
            -12, 12, 0
            -7.8, 2.4, 6.64
            -7.8, 7.2,  6
            -7.8, 12, 2.8
            -3.6, 2.4, 7.48
            -3.6, 7.2, 6.84
            -3.6, 12, 3.64    ];
        % weights
        W = [
            1
            1
            1
            1
            1
            1
            1
            1
            1
            ];
        % knots
        U = [
            0
            0
            1
            1
            ];
        
        V = [
            0
            0
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(4,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
    case 4  % torus 4
        Npat = 4;
        multiPatSurf = cell(Npat,7); % p,q,m,n,U,V,Ctrlpts.
        
        %% patch 1
        p = 2;  q = 2;
        m = 4;  n = 4;
        % contorl point
        Points = [
            70, 0, 0
            70, 0, 20
            50, 0, 20
            30, 0, 20
            30, 0, 0
            70, 70, 0
            70, 70, 20
            50, 50, 20
            30, 30, 20
            30, 30, 0
            0, 70, 0
            0, 70, 20
            0, 50, 20
            0, 30, 20
            0, 30, 0
            -70, 70, 0
            -70, 70, 20
            -50, 50, 20
            -30, 30, 20
            -30, 30, 0
            -70, 0, 0
            -70, 0, 20
            -50, 0, 20
            -30, 0, 20
            -30, 0, 0 ];
        % weights
        W = [
            1
            0.707107
            1
            0.707107
            1
            0.707107
            0.5
            0.707107
            0.5
            0.707107
            1
            0.707107
            1
            0.707107
            1
            0.707107
            0.5
            0.707107
            0.5
            0.707107
            1
            0.707107
            1
            0.707107
            1
            ];
        % knots
        U = [
            0
            0
            0.5
            0.5
            1
            1
            ];
        
        V = [
            0
            0
            0.5
            0.5
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(1,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 2
        p = 2;  q = 2;
        m = 4;  n = 4;
        % contorl point
        Points = [
            30, 0, 0
            30, 0, -20
            50, 0, -20
            70, 0, -20
            70, 0, 0
            30, 30, 0
            30, 30, -20
            50, 50, -20
            70, 70, -20
            70, 70, 0
            0, 30, 0
            0, 30, -20
            0, 50, -20
            0, 70, -20
            0, 70, 0
            -30, 30, 0
            -30, 30, -20
            -50, 50, -20
            -70, 70, -20
            -70, 70, 0
            -30, 0, 0
            -30, 0, -20
            -50, 0, -20
            -70, 0, -20
            -70, 0, 0  ];
        % weights
        W = [
            1
            0.707107
            1
            0.707107
            1
            0.707107
            0.5
            0.707107
            0.5
            0.707107
            1
            0.707107
            1
            0.707107
            1
            0.707107
            0.5
            0.707107
            0.5
            0.707107
            1
            0.707107
            1
            0.707107
            1
            ];
        % knots
        U = [
            0
            0
            0.5
            0.5
            1
            1
            ];
        
        V = [
            0
            0
            0.5
            0.5
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(2,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 3
        p = 2;  q = 2;
        m = 4;  n = 4;
        % contorl point
        Points = [
            -30, 0, 0
            -30, 0, -20
            -50, 0, -20
            -70, 0, -20
            -70, 0, 0
            -30, -30, 0
            -30, -30, -20
            -50, -50, -20
            -70, -70, -20
            -70, -70, 0
            0, -30, 0
            0, -30, -20
            0, -50, -20
            0, -70, -20
            0, -70, 0
            30, -30, 0
            30, -30, -20
            50, -50, -20
            70, -70, -20
            70, -70, 0
            30, 0, 0
            30, 0, -20
            50, 0, -20
            70, 0, -20
            70, 0, 0  ];
        % weights
        W = [
            1
            0.707107
            1
            0.707107
            1
            0.707107
            0.5
            0.707107
            0.5
            0.707107
            1
            0.707107
            1
            0.707107
            1
            0.707107
            0.5
            0.707107
            0.5
            0.707107
            1
            0.707107
            1
            0.707107
            1
            ];
        % knots
        U = [
            0
            0
            0.5
            0.5
            1
            1
            ];
        
        V = [
            0
            0
            0.5
            0.5
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(3,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 4
        p = 2;  q = 2;
        m = 4;  n = 4;
        % contorl point
        Points = [
            -70, 0, 0
            -70, 0, 20
            -50, 0, 20
            -30, 0, 20
            -30, 0, 0
            -70, -70, 0
            -70, -70, 20
            -50, -50, 20
            -30, -30, 20
            -30, -30, 0
            0, -70, 0
            0, -70, 20
            0, -50, 20
            0, -30, 20
            0, -30, 0
            70, -70, 0
            70, -70, 20
            50, -50, 20
            30, -30, 20
            30, -30, 0
            70, 0, 0
            70, 0, 20
            50, 0, 20
            30, 0, 20
            30, 0, 0  ];
        % weights
        W = [
            1
            0.707107
            1
            0.707107
            1
            0.707107
            0.5
            0.707107
            0.5
            0.707107
            1
            0.707107
            1
            0.707107
            1
            0.707107
            0.5
            0.707107
            0.5
            0.707107
            1
            0.707107
            1
            0.707107
            1
            ];
        % knots
        U = [
            0
            0
            0.5
            0.5
            1
            1
            ];
        
        V = [
            0
            0
            0.5
            0.5
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(4,:) = {p,q,m,n,U,V,Ctrlpts} ;
                
    case 5  % tablet holder -- 4 patches
        Npat = 4;
        multiPatSurf = cell(Npat,7); % p,q,m,n,U,V,Ctrlpts.
        
        %% patch 1
        p = 3;  q = 3;
        m = 4;  n = 4;
        % contorl point
        Points = [
            0, 0, 6
            0, 1.0, 6.0
            0, 3.0, 4.0
            0, 6.0, 1.0
            0, 9, 0
            1, 0, 6
            1, 1.0, 5.666667
            1, 3.0, 3.666667
            1.0, 6.0, 1.0
            1.0, 9, 0
            3, 0, 4
            3.0, 1.0, 3.666667
            3.0, 3.0, 3.0
            3.0, 6.0, 2.333333
            3, 9, 2
            6, 0, 1
            6.0, 1.0, 1.0
            6, 3.0, 2.333333
            6, 6.0, 4.333333
            6, 9, 5
            9, 0, 0
            9, 1, 0
            9, 3, 2
            9, 6, 5
            9, 9, 6  ];
        % weights
        W = [
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1 ];
        % knots
        U = [
            0
            0
            0
            0.5
            1
            1
            1
            ];
        
        V = [
            0
            0
            0
            0.5
            1
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(1,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 2
        p = 3;  q = 3;
        m = 4;  n = 4;
        % contorl point
        Points = [
            0, 0, 6
            -1, 0, 6
            -3, 0, 4
            -6, 0, 1
            -9, 0, 0
            0, 1.0, 6.0
            -1.0, 1.0, 5.666667
            -3.0, 1.0, 3.666667
            -6.0, 1.0, 1.0
            -9, 1, 0
            0, 3.0, 4.0
            -1, 3.0, 3.666667
            -3.0, 3.0, 3.0
            -6.0, 3.0, 2.333333
            -9, 3, 2
            0, 6.0, 1
            -1, 6.0, 1.0
            -3, 6.0, 2.333333
            -6, 6.0, 4.333333
            -9, 6, 5
            0, 9, 0
            -1.0, 9, 0
            -3, 9, 2
            -6, 9, 5
            -9, 9, 6 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1 ];
        % knots
        U = [
            0
            0
            0
            0.5
            1
            1
            1
            ];
        
        V = [
            0
            0
            0
            0.5
            1
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(2,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 3
        p = 3;  q = 3;
        m = 4;  n = 4;
        % contorl point
        Points = [
            0, 0, 6
            0, -1.0, 6.0
            0, -3.0, 4.0
            0, -6.0, 1
            0, -9, 0
            -1, 0, 6
            -1.0, -1.0, 5.666667
            -1.0, -3.0, 3.666667
            -1, -6.0, 1.0
            -1.0, -9, 0
            -3, 0, 4
            -3.0, -1.0, 3.666667
            -3, -3.0, 3.0
            -3, -6.0, 2.333333
            -3, -9, 2
            -6, 0, 1
            -6.0, -1.0, 1.0
            -6, -3.0, 2.333333
            -6, -6.0, 4.333333
            -6, -9, 5
            -9, 0, 0
            -9, -1, 0
            -9, -3, 2
            -9, -6, 5
            -9, -9, 6   ];
        % weights
        W = [
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1 ];
        % knots
        U = [
            0
            0
            0
            0.5
            1
            1
            1
            ];
        
        V = [
            0
            0
            0
            0.5
            1
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(3,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 4
        p = 3;  q = 3;
        m = 4;  n = 4;
        % contorl point
        Points = [
            0, 0, 6
            1, 0, 6
            3, 0, 4
            6, 0, 1
            9, 0, 0
            0, -1.0, 6.0
            1.0, -1.0, 5.666667
            3, -1.0, 3.666667
            6, -1.0, 1.0
            9, -1, 0
            0, -3.0, 4.0
            1, -3.0, 3.666667
            3, -3.0, 3.0
            6, -3.0, 2.333333
            9, -3, 2
            0, -6.0, 1
            1, -6.0, 1.0
            3, -6.0, 2.333333
            6, -6.0, 4.333333
            9, -6, 5
            0, -9, 0
            1.0, -9, 0
            3, -9, 2
            6, -9, 5
            9, -9, 6 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1
            1 ];
        % knots
        U = [
            0
            0
            0
            0.5
            1
            1
            1
            ];
        
        V = [
            0
            0
            0
            0.5
            1
            1
            1
            ];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(4,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
    case 6  % pipes -- 24 patches
        Npat = 24;
        multiPatSurf = cell(Npat,7); % p,q,m,n,U,V,Ctrlpts.
        
        %% patch 1
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            0, -10, 5
            0, -24, 5
            5, -10, 5
            5, -24, 5
            5, -10, 0
            5, -24, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(1,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 2
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            0, -10, -5
            0, -24, -5
            5, -10, -5
            5, -24, -5
            5, -10, 0
            5, -24, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(2,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 3
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            0, -10, -5
            0, -24, -5
            -5, -10, -5
            -5, -24, -5
            -5, -10, 0
            -5, -24, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(3,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 4
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            0, -10, 5
            0, -24, 5
            -5, -10, 5
            -5, -24, 5
            -5, -10, 0
            -5, -24, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(4,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 5
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, 5
            0, -3.333333, 5
            0, -6.666667, 5
            0, -10, 5
            3.861696, -2.229551, 5
            3.411035, -4.385139, 5.0
            3.074878, -6.96837, 5.0
            2.928932, -10.0, 5.0
            6.655767, -3.842709, 3.936742
            5.81463, -5.461531, 3.647581
            5.211507, -7.532714, 3.315144
            5, -10, 2.928932
            6.655767, -3.842709, 0
            5.59049, -5.687822, 0
            5.0, -7.807054, 0
            5, -10, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(5,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 6
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, -5
            0, -3.333333, -5
            0, -6.666667, -5
            0, -10, -5
            3.861696, -2.229551, -5
            3.411035, -4.385139, -5.0
            3.074878, -6.96837, -5.0
            2.928932, -10.0, -5.0
            6.655767, -3.842709, -3.936742
            5.81463, -5.461531, -3.647581
            5.211507, -7.532714, -3.315144
            5, -10, -2.928932
            6.655767, -3.842709, 0
            5.59049, -5.687822, 0
            5.0, -7.807054, 0
            5, -10, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(6,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 7
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, -5
            0, -3.333333, -5
            0, -6.666667, -5
            0, -10, -5
            -3.861696, -2.229551, -5
            -3.411035, -4.385139, -5.0
            -3.074878, -6.96837, -5.0
            -2.928932, -10.0, -5.0
            -6.655767, -3.842709, -3.936742
            -5.81463, -5.461531, -3.647581
            -5.211507, -7.532714, -3.315144
            -5, -10, -2.928932
            -6.655767, -3.842709, 0
            -5.59049, -5.687822, 0
            -5.0, -7.807054, 0
            -5, -10, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(7,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 8
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, 5
            0, -3.333333, 5
            0, -6.666667, 5
            0, -10, 5
            -3.861696, -2.229551, 5
            -3.411035, -4.385139, 5.0
            -3.074878, -6.96837, 5.0
            -2.928932, -10.0, 5.0
            -6.655767, -3.842709, 3.936742
            -5.81463, -5.461531, 3.647581
            -5.211507, -7.532714, 3.315144
            -5, -10, 2.928932
            -6.655767, -3.842709, 0
            -5.59049, -5.687822, 0
            -5.0, -7.807054, 0
            -5, -10, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(8,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        
        %% patch 9
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            8.660254, 5.0, 5
            20.78461, 12.0, 5
            6.160254, 9.330127, 5
            18.28461, 16.330127, 5
            6.160254, 9.330127, 0
            18.28461, 16.330127, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(9,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 10
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            8.660254, 5.0, -5
            20.78461, 12.0, -5
            6.160254, 9.330127, -5
            18.28461, 16.330127, -5
            6.160254, 9.330127, 0
            18.28461, 16.330127, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(10,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 11
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            8.660254, 5.0, -5
            20.78461, 12.0, -5
            11.160254, 0.669873, -5
            23.28461, 7.669873, -5
            11.160254, 0.669873, 0
            23.28461, 7.669873, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(11,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 12
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            8.660254, 5.0, 5
            20.78461, 12.0, 5
            11.160254, 0.669873, 5
            23.28461, 7.669873, 5
            11.160254, 0.669873, 0
            23.28461, 7.669873, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(12,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 13
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, 5
            2.886751, 1.666667, 5
            5.773503, 3.333333, 5
            8.660254, 5.0, 5
            0, 4.459103, 5
            2.092124, 5.146613, 5.0
            4.497347, 6.147107, 5.0
            7.195788, 7.53653, 5.0
            0, 7.685417, 3.936742
            1.822509, 7.766383, 3.647581
            3.917768, 8.279654, 3.315144
            6.160254, 9.330127, 2.928932
            0, 7.685417, 0
            2.130553, 7.685417, 0
            4.261107, 8.233654, 0
            6.160254, 9.330127, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(13,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 14
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, -5
            2.886751, 1.666667, -5
            5.773503, 3.333333, -5
            8.660254, 5.0, -5
            0, 4.459103, -5
            2.092124, 5.146613, -5.0
            4.497347, 6.147107, -5.0
            7.195788, 7.53653, -5.0
            0, 7.685417, -3.936742
            1.822509, 7.766383, -3.647581
            3.917768, 8.279654, -3.315144
            6.160254, 9.330127, -2.928932
            0, 7.685417, 0
            2.130553, 7.685417, 0
            4.261107, 8.233654, 0
            6.160254, 9.330127, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(14,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 15
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, -5
            2.886751, 1.666667, -5
            5.773503, 3.333333, -5
            8.660254, 5.0, -5
            3.861696, -2.229551, -5
            5.503159, -0.761474, -5.0
            7.572225, 0.821263, -5.0
            10.12472, 2.46347, -5.0
            6.655767, -3.842709, -3.936742
            7.63714, -2.304852, -3.647581
            9.129275, -0.74694, -3.315144
            11.160254, 0.669873, -2.928932
            6.655767, -3.842709, 0
            7.721043, -1.997595, 0
            9.261107, -0.4266, 0
            11.160254, 0.669873, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(15,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 16
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, 5
            2.886751, 1.666667, 5
            5.773503, 3.333333, 5
            8.660254, 5.0, 5
            3.861696, -2.229551, 5
            5.503159, -0.761474, 5.0
            7.572225, 0.821263, 5.0
            10.12472, 2.46347, 5.0
            6.655767, -3.842709, 3.936742
            7.63714, -2.304852, 3.647581
            9.129275, -0.74694, 3.315144
            11.160254, 0.669873, 2.928932
            6.655767, -3.842709, 0
            7.721043, -1.997595, 0
            9.261107, -0.4266, 0
            11.160254, 0.669873, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(16,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        
        %% patch 17
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            -8.660254, 5, 5
            -20.78461, 12, 5
            -11.160254, 0.669873, 5
            -23.28461, 7.669873, 5
            -11.160254, 0.669873, 0
            -23.28461, 7.669873, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(17,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 18
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            -8.660254, 5, -5
            -20.78461, 12, -5
            -11.160254, 0.669873, -5
            -23.28461, 7.669873, -5
            -11.160254, 0.669873, 0
            -23.28461, 7.669873, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(18,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 19
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            -8.660254, 5, -5
            -20.78461, 12, -5
            -6.160254, 9.330127, -5
            -18.28461, 16.330127, -5
            -6.160254, 9.330127, 0
            -18.28461, 16.330127, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(19,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 20
        p = 2;  q = 1;
        m = 2;  n = 1;
        % contorl point
        Points = [
            -8.660254, 5, 5
            -20.78461, 12, 5
            -6.160254, 9.330127, 5
            -18.28461, 16.330127, 5
            -6.160254, 9.330127, 0
            -18.28461, 16.330127, 0 ];
        % weights
        W = [
            1
            1
            0.707107
            0.707107
            1
            1];
        % knots
        U = [
            0
            0
            1
            1];
        
        V = [
            0
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(20,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 21
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, 5
            -2.886751, 1.666667, 5
            -5.773503, 3.333333, 5
            -8.660254, 5, 5
            -3.861696, -2.229551, 5
            -5.503159, -0.761474, 5.0
            -7.572225, 0.821263, 5.0
            -10.12472, 2.46347, 5.0
            -6.655767, -3.842709, 3.936742
            -7.63714, -2.304852, 3.647581
            -9.129275, -0.74694, 3.315144
            -11.160254, 0.669873, 2.928932
            -6.655767, -3.842709, 0
            -7.721043, -1.997595, 0
            -9.261107, -0.4266, 0
            -11.160254, 0.669873, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(21,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 22
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, -5
            -2.886751, 1.666667, -5
            -5.773503, 3.333333, -5
            -8.660254, 5, -5
            -3.861696, -2.229551, -5
            -5.503159, -0.761474, -5.0
            -7.572225, 0.821263, -5.0
            -10.12472, 2.46347, -5.0
            -6.655767, -3.842709, -3.936742
            -7.63714, -2.304852, -3.647581
            -9.129275, -0.74694, -3.315144
            -11.160254, 0.669873, -2.928932
            -6.655767, -3.842709, 0
            -7.721043, -1.997595, 0
            -9.261107, -0.4266, 0
            -11.160254, 0.669873, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(22,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 23
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, -5
            -2.886751, 1.666667, -5
            -5.773503, 3.333333, -5
            -8.660254, 5, -5
            0, 4.459103, -5
            -2.092124, 5.146613, -5.0
            -4.497347, 6.147107, -5.0
            -7.195788, 7.53653, -5.0
            0, 7.685417, -3.936742
            -1.822509, 7.766383, -3.647581
            -3.917768, 8.279654, -3.315144
            -6.160254, 9.330127, -2.928932
            0, 7.685417, 0
            -2.130553, 7.685417, 0
            -4.261107, 8.233654, 0
            -6.160254, 9.330127, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(23,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 24
        p = 3;  q = 3;
        m = 3;  n = 3;
        % contorl point
        Points = [
            0, 0, 5
            -2.886751, 1.666667, 5
            -5.773503, 3.333333, 5
            -8.660254, 5, 5
            0, 4.459103, 5
            -2.092124, 5.146613, 5.0
            -4.497347, 6.147107, 5.0
            -7.195788, 7.53653, 5.0
            0, 7.685417, 3.936742
            -1.822509, 7.766383, 3.647581
            -3.917768, 8.279654, 3.315144
            -6.160254, 9.330127, 2.928932
            0, 7.685417, 0
            -2.130553, 7.685417, 0
            -4.261107, 8.233654, 0
            -6.160254, 9.330127, 0 ];
        % weights
        W = [
            1
            1
            1
            1
            1
            0.934913
            0.869825
            0.804738
            1
            0.934913
            0.869825
            0.804738
            1
            1
            1
            1];
        % knots
        U = [
            0
            0
            0
            1
            1
            1];
        
        V = [
            0
            0
            0
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(24,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
    case 3  % bird's nest -- 4 patches
        Npat = 4;
        multiPatSurf = cell(Npat,7); % p,q,m,n,U,V,Ctrlpts.
        
        %% patch 1
        p = 3;  q = 3;
        m = 4;  n = 7;
        % contorl point
        Points = [
            45, 0, 0
            45, 2.943373, 0
            43.664387, 9.075392, 0
            37.24989, 17.787705, 0
            26.681557, 24.83326, 0
            13.613089, 29.109591, 0
            4.415059, 30.0, 0
            0, 30, 0
            55.566794, 0, 23.987119
            55.517147, 4.612307, 24.399491
            52.601956, 13.918603, 24.542752
            42.57371, 24.475431, 23.335096
            29.687411, 31.732698, 21.361694
            15.295387, 35.966885, 18.976258
            5.273717, 36.940721, 17.29756
            0.277398, 36.982823, 16.483376
            60.802785, 0, 45.014822
            60.442123, 5.598623, 44.589895
            56.187431, 16.51125, 43.079438
            44.098205, 27.603643, 39.918355
            30.302116, 34.703107, 36.813933
            15.696583, 38.946451, 34.089681
            5.586395, 40.216327, 32.708086
            0.411383, 40.471738, 32.178018
            43.774562, 0, 32.146618
            44.005577, 3.593148, 32.256266
            42.207114, 10.981264, 32.05612
            34.569563, 19.650309, 30.93306
            24.107653, 25.647698, 29.487017
            12.203272, 28.94009, 28.02587
            3.981262, 29.443905, 27.193184
            -0.045428, 29.300358, 26.84687
            33.0, 0, 21
            33.0, 2.180276, 21
            32.010657, 6.722513, 21.0
            27.259177, 13.176078, 21.0
            19.430783, 18.395007, 21
            9.750436, 21.56266, 21
            2.937081, 22.222222, 21
            -0.333333, 22.222222, 21 ];
        % weights
        W = [
            1
            0.960948
            0.898464
            0.851601
            0.851601
            0.898464
            0.960948
            1
            1.0
            0.992396
            0.98023
            0.971105
            0.971105
            0.98023
            0.992396
            1.0
            1
            1.014747
            1.038343
            1.056039
            1.056039
            1.038343
            1.014747
            1
            1
            0.983159
            0.956214
            0.936004
            0.936004
            0.956214
            0.983159
            1
            1
            0.960948
            0.898464
            0.851601
            0.851601
            0.898464
            0.960948
            1];
        % knots
        U = [
            0
            0
            0
            0.589689
            1
            1
            1];
        
        V = [
            0
            0
            0
            0.2
            0.4
            0.6
            0.8
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(1,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 2
        p = 3;  q = 3;
        m = 4;  n = 7;
        % contorl point
        Points = [
            0, 30, 0
            -4.415059, 30, 0
            -13.613089, 29.109591, 0
            -26.681557, 24.83326, 0
            -37.24989, 17.787705, 0
            -43.664387, 9.075392, 0
            -45.0, 2.943373, 0
            -45, 0, 0
            0.277398, 36.982823, 16.483376
            -4.721926, 36.962615, 16.378945
            -14.793536, 36.033379, 16.186207
            -29.360187, 31.8121, 16.225291
            -42.400903, 24.496398, 16.844726
            -52.438063, 13.910556, 17.863721
            -55.338592, 4.609655, 17.907791
            -55.382392, -2.5286e-7, 17.544782
            0.411383, 40.471738, 32.178018
            -4.782864, 40.252955, 31.171206
            -14.995422, 39.053839, 29.583815
            -29.889227, 34.828017, 28.733694
            -43.928235, 27.636626, 29.708146
            -56.032065, 16.498254, 32.292961
            -60.26303, 5.594186, 33.729026
            -60.61366, -4.3257e-7, 33.993786
            -0.045428, 29.300358, 26.84687
            -4.059673, 29.459513, 26.538288
            -12.293467, 28.988233, 26.005818
            -24.307045, 25.705881, 25.723233
            -34.882107, 19.665673, 26.177134
            -42.541998, 10.975438, 27.220371
            -44.350074, 3.591258, 27.628227
            -44.12702, -1.7859e-7, 27.59652
            -0.333333, 22.222222, 21
            -3.603748, 22.222222, 21
            -10.417103, 21.56266, 21
            -20.09745, 18.395007, 21
            -27.925844, 13.176078, 21
            -32.677324, 6.722513, 21
            -33.666667, 2.180276, 21.0
            -33.666667, 0, 21  ];
        % weights
        W = [
            1
            0.960948
            0.898464
            0.851601
            0.851601
            0.898464
            0.960948
            1
            1.0
            0.992396
            0.98023
            0.971105
            0.971105
            0.98023
            0.992396
            1.0
            1
            1.014747
            1.038343
            1.056039
            1.056039
            1.038343
            1.014747
            1
            1
            0.983159
            0.956214
            0.936004
            0.936004
            0.956214
            0.983159
            1
            1
            0.960948
            0.898464
            0.851601
            0.851601
            0.898464
            0.960948
            1];
        % knots
        U = [
            0
            0
            0
            0.589689
            1
            1
            1];
        
        V = [
            0
            0
            0
            0.2
            0.4
            0.6
            0.8
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(2,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 3
        p = 3;  q = 3;
        m = 4;  n = 7;
        % contorl point
        Points = [
            -45, 0, 0
            -45.0, -2.943373, 0
            -43.664387, -9.075392, 0
            -37.24989, -17.787705, 0
            -26.681557, -24.83326, 0
            -13.613089, -29.109591, 0
            -4.415059, -30, 0
            0, -30, 0
            -55.382392, -2.5286e-7, 17.544782
            -55.338577, -4.609656, 17.907791
            -52.438019, -13.910528, 17.863721
            -42.400909, -24.496379, 16.844613
            -29.360203, -31.812148, 16.225427
            -14.79352, -36.033354, 16.186098
            -4.721923, -36.962582, 16.378914
            0.277399, -36.9828, 16.483365
            -60.61366, -4.3257e-7, 33.993786
            -60.263006, -5.594187, 33.729025
            -56.031993, -16.498209, 32.292959
            -43.928244, -27.636597, 29.707969
            -29.889253, -34.828093, 28.733908
            -14.995397, -39.053798, 29.583638
            -4.78286, -40.252901, 31.171155
            0.411385, -40.471699, 32.178
            -44.12702, -1.7859e-7, 27.59652
            -44.350064, -3.591258, 27.628227
            -42.541966, -10.975418, 27.220371
            -34.882111, -19.665659, 26.177051
            -24.307057, -25.705916, 25.723332
            -12.293456, -28.988215, 26.005739
            -4.059671, -29.45949, 26.538266
            -0.045427, -29.300342, 26.846862
            -33.666667, 0, 21
            -33.666667, -2.180276, 21.0
            -32.677324, -6.722513, 21.0
            -27.925844, -13.176078, 21.0
            -20.09745, -18.395007, 21
            -10.417103, -21.56266, 21
            -3.603748, -22.222222, 21
            -0.333333, -22.222222, 21  ];
        % weights
        W = [
            1
            0.960948
            0.898464
            0.851601
            0.851601
            0.898464
            0.960948
            1
            1.0
            0.992396
            0.98023
            0.971105
            0.971105
            0.98023
            0.992396
            1.0
            1
            1.014747
            1.038343
            1.056039
            1.056039
            1.038343
            1.014747
            1
            1
            0.983159
            0.956214
            0.936004
            0.936004
            0.956214
            0.983159
            1
            1
            0.960948
            0.898464
            0.851601
            0.851601
            0.898464
            0.960948
            1];
        % knots
        U = [
            0
            0
            0
            0.589689
            1
            1
            1];
        
        V = [
            0
            0
            0
            0.2
            0.4
            0.6
            0.8
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(3,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
        %% patch 4
        p = 3;  q = 3;
        m = 4;  n = 7;
        % contorl point
        Points = [
            0, -30, 0
            4.415059, -30, 0
            13.613089, -29.109591, 0
            26.681557, -24.83326, 0
            37.24989, -17.787705, 0
            43.664387, -9.075392, 0
            45.0, -2.943373, 0
            45, 0, 0
            0.277399, -36.9828, 16.483365
            5.273717, -36.940708, 17.297568
            15.295393, -35.966919, 18.976266
            29.687414, -31.732706, 21.361689
            42.573718, -24.475435, 23.335099
            52.601932, -13.91859, 24.542748
            55.51714, -4.612307, 24.399489
            55.566794, 0, 23.987119
            0.411385, -40.471699, 32.178
            5.586394, -40.216305, 32.7081
            15.696593, -38.946506, 34.089694
            30.302121, -34.70312, 36.813925
            44.098217, -27.603649, 39.91836
            56.187392, -16.511229, 43.079431
            60.44211, -5.598623, 44.589893
            60.802785, 0, 45.014822
            -0.045427, -29.300342, 26.846862
            3.981261, -29.443896, 27.19319
            12.203276, -28.940114, 28.025876
            24.107656, -25.647704, 29.487013
            34.569569, -19.650312, 30.933063
            42.207097, -10.981255, 32.056118
            44.005571, -3.593148, 32.256265
            43.774562, 0, 32.146618
            -0.333333, -22.222222, 21
            2.937081, -22.222222, 21
            9.750436, -21.56266, 21
            19.430783, -18.395007, 21
            27.259177, -13.176078, 21
            32.010657, -6.722513, 21
            33.0, -2.180276, 21.0
            33.0, 0, 21  ];
        % weights
        W = [
            1
            0.960948
            0.898464
            0.851601
            0.851601
            0.898464
            0.960948
            1
            1.0
            0.992396
            0.98023
            0.971105
            0.971105
            0.98023
            0.992396
            1.0
            1
            1.014747
            1.038343
            1.056039
            1.056039
            1.038343
            1.014747
            1
            1
            0.983159
            0.956214
            0.936004
            0.936004
            0.956214
            0.983159
            1
            1
            0.960948
            0.898464
            0.851601
            0.851601
            0.898464
            0.960948
            1];
        % knots
        U = [
            0
            0
            0
            0.589689
            1
            1
            1];
        
        V = [
            0
            0
            0
            0.2
            0.4
            0.6
            0.8
            1
            1
            1];
        % save
        U = [U(1); U; U(end)]';
        V = [V(1); V; V(end)]';
        Ctrlpts = zeros(m+1,n+1,3);
        Weights = ones(m+1,n+1);
        
        for i = 1:m+1
            for j = 1:n+1
                Ctrlpts(i,j,:) = Points((i-1)*(n+1)+j,:);
                Weights(i,j) = W((i-1)*(n+1)+j,:);
            end
        end
        
        Ctrlpts(:,:,4) = Weights;
        multiPatSurf(4,:) = {p,q,m,n,U,V,Ctrlpts} ;
        
end

end