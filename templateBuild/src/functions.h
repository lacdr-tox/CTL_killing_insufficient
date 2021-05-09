//
// Created by richard on 10/01/2019.
//
// return ratio of ctl to target cells
float getETRatio(float time, int method){
    float ETRatio;
    if(method ==0){// from Bousso paper
        if(time<5.) return(0);
        if(time>=5. & time<7.) return (time-5.)*0.5*12500./1000000.;
        if(time>=7.) return ((time-7.)*12500.+12500.)/1000000.;
    }
    if(method==1){
        float insertOnDay = 5.;
        float peak=2.5;
        float rIncrease=2.5;
        float initRatio = 0.0001;
        if(time>insertOnDay+peak) return(0);
        // since this is infrequently called it is k to get the exp here:
        ETRatio = initRatio*exp(rIncrease*(time-insertOnDay));

    }
    return(ETRatio);
}

