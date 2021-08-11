


for i=1:Ttime
    GLCe_current = GLCe_next;
    GLCn_current = GLCn_next;
    
    v_en = Ven_max * ((GLCe_current/(GLCe_current + Ken_T)) - (GLCn_current/(GLCn_current + Kce_T)));
    
    v_ce = Vce_max * ((GLCc_current/(GLCc_current + Kce_T)) - (GLCe_current/(GLCe_current + Kce_T)));

    fG6Pn = 1/(1 + exp(-20*(t-.6)));  % .6 is an event time, but what event, TODO. Note: SS value of fG6Pn = 0.75
    vn_hk = kn_hk * ATP * (GLCn_current/(GLCn_current + Km))*(1 - fG6Pn); 
    
    changeGLCe = v_ce - R_ne*v_en; % ignoring astro

    changeGCLn = v_en - vn_hk;  % why R_ne*v_en leaving extra fluid for neuron but v_en entering neuron from fluid TODO
    
    GLCe_next = GLCe_current + changeGLCe*dt; % TODO: the GLCe variable will be used in the equation for v_hk in the Muddapu code
    GLCn_next = GLCn_current + changeGCLn*dt;
    
    
    
end