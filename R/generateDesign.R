ones<-function(m,n)
  b=matrix(rep(1,m*n),nrow=m)



generateDesign<-function(pattern,pr)
{
  beta1=beta2=beta3=matrix(rep(0,64*64),nrow=64);
  prmter = pr*c(0,1,2,3,4)/5;
  if (pattern == 1)
  {
    pattern1 = matrix(rep(1,24*24),nrow=24)*prmter[1];
    pattern1[1:12,1:12] = matrix(rep(1,12*12),nrow=12)*prmter[3];
    pattern1[1:12,13:24] = matrix(rep(1,12*12),nrow=12)*prmter[5];
    pattern1[13:24,1:12] = matrix(rep(1,12*12),nrow=12)*prmter[4];
    pattern1[13:24,13:24] = matrix(rep(1,12*12),nrow=12)*prmter[2];
    
    #pattern 2 %%%%
      pattern2 = matrix(rep(1,24*24),nrow=24)*prmter[1];
      pattern2[1:16,1:8] = matrix(rep(1,16*8),nrow=16)*prmter[2];
      pattern2[1:8,9:24] = matrix(rep(1,16*8),nrow=8)*prmter[3];
      pattern2[17:24,1:16] = matrix(rep(1,16*8),nrow=8)*prmter[4];
      pattern2[9:24,17:24] = matrix(rep(1,16*8),nrow=16)*prmter[5];
      
      #pattern 3 %%%%
        pattern3 = matrix(rep(1,24*24),nrow=24)*prmter[1];
        pattern3[1:24,1:6] = matrix(rep(1,24*6),nrow=24)*prmter[3];
        pattern3[1:24,7:12] = matrix(rep(1,24*6),nrow=24)*prmter[5];
        pattern3[1:24,13:18] = matrix(rep(1,24*6),nrow=24)*prmter[2];
        pattern3[1:24,19:24] = matrix(rep(1,24*6),nrow=24)*prmter[4];
        
        
       # pattern 4 %%%%
          pattern4 = matrix(rep(1,24*24),nrow=24)*prmter[1];
          pattern4[1:24,1:4] = matrix(rep(1,24*4),nrow=24)*prmter[5];
          pattern4[1:24,5:8] = matrix(rep(1,24*4),nrow=24)*prmter[3];
          pattern4[1:24,9:12] = matrix(rep(1,24*4),nrow=24)*prmter[4];
          pattern4[1:24,13:16] = matrix(rep(1,24*4),nrow=24)*prmter[3];
          pattern4[1:24,17:20] = matrix(rep(1,24*4),nrow=24)*prmter[2];
          pattern4[1:24,21:24] = matrix(rep(1,24*4),nrow=24)*prmter[4];
          
          # over all pattern %%%%
            beta1[5:28,5:28] = t(pattern4);
            beta1[5:28,37:60] = pattern3;
            beta1[37:60,5:28] = pattern2;
            beta1[37:60,37:60] = t(pattern1);
            
            beta2[5:28,5:28] = pattern1;
            beta2[5:28,37:60] = pattern2;
            beta2[37:60,5:28] = pattern3;
            beta2[37:60,37:60] = t(pattern4);
            
            beta3[5:28,5:28] = pattern3;
            beta3[5:28,37:60] = t(pattern4);
            beta3[37:60,5:28] = pattern1;
            beta3[37:60,37:60] = t(pattern2);
    
    
  }
  
  if(pattern == 2)
  {
    
    
      pattern1 = matrix(rep(1,18*18),nrow=18)*prmter[1];
      pattern1[1:9,1:9] = matrix(rep(1,81),nrow=9)*prmter[3];
      pattern1[1:9,10:18] = matrix(rep(1,81),nrow=9)*prmter[5];
      pattern1[10:18,1:9] = matrix(rep(1,81),nrow=9)*prmter[4];
      pattern1[10:18,10:18] = matrix(rep(1,81),nrow=9)*prmter[2];
      
  
        pattern2 = matrix(rep(1,18*18),nrow=18)*prmter[1];
        pattern2[1:12,1:6] = matrix(rep(1,72),nrow=12)*prmter[2];
        pattern2[1:6,7:18] = matrix(rep(1,72),nrow=6)*prmter[3];
        pattern2[13:18,1:12] = matrix(rep(1,72),nrow=6)*prmter[4];
        pattern2[7:18,13:18] = matrix(rep(1,72),nrow=12)*prmter[5];
        
        
          pattern3 = matrix(rep(1,18*18),nrow=18)*prmter[1];
          pattern3[1:18,1:6] = matrix(rep(1,108),nrow=18)*prmter[3];
          pattern3[1:18,7:12] = matrix(rep(1,108),nrow=18)*prmter[4];
          pattern3[1:18,13:18] = matrix(rep(1,108),nrow=18)*prmter[2];
          
         
            pattern4 = matrix(rep(1,18*18),nrow=18)*prmter[1];
            pattern4[1:18,1:3] = matrix(rep(1,54),nrow=18)*prmter[2];
            pattern4[1:18,4:6] = matrix(rep(1,54),nrow=18)*prmter[4];
            pattern4[1:18,7:9] = matrix(rep(1,54),nrow=18)*prmter[5];
            pattern4[1:18,10:12] = matrix(rep(1,54),nrow=18)*prmter[3];
            pattern4[1:18,13:15] = matrix(rep(1,54),nrow=18)*prmter[4];
            pattern4[1:18,16:18] = matrix(rep(1,54),nrow=18)*prmter[2];
            
              beta1[9:26,9:26] = pattern1;
              beta1[9:26,41:58] = pattern3;
              beta1[41:58,9:26] = pattern2;
              beta1[41:58,41:58] = t(pattern4);
              
              beta2[9:26,9:26] = pattern1;
              beta2[9:26,41:58] = pattern2;
              beta2[41:58,9:26] = pattern3;
              beta2[41:58,41:58] = t(pattern4);
              
              beta3[9:26,9:26] = pattern3;
              beta3[9:26,41:58] = t(pattern4);
              beta3[41:58,9:26] = pattern1;
              beta3[41:58,41:58] = t(pattern2);
  }
if(pattern == 3)
{
  pattern1 = ones(16,16)*prmter[1];
  pattern1[1:8,1:8] = ones(8,8)*prmter[3];
  pattern1[1:8,9:16] = ones(8,8)*prmter[5];
  pattern1[9:16,1:8] = ones(8,8)*prmter[4];
  pattern1[9:16,9:16] = ones(8,8)*prmter[2];
  
    pattern2 = ones(16,16)*prmter[1];
    pattern2[1:10,1:6] = ones(10,6)*prmter[2];
    pattern2[1:6,7:16] = ones(6,10)*prmter[3];
    pattern2[11:16,1:10] = ones(6,10)*prmter[4];
    pattern2[7:16,11:16] = ones(10,6)*prmter[5];
    
      pattern3 = ones(16,16)*prmter[1];
      pattern3[1:16,1:4] = ones(16,4)*prmter[3];
      pattern3[1:16,5:8] = ones(16,4)*prmter[5];
      pattern3[1:16,9:12] = ones(16,4)*prmter[4];
      pattern3[1:16,13:16] = ones(16,4)*prmter[2];
      
        pattern4 = ones(16,16)*prmter[1];
        pattern4[1:16,1:2] = ones(16,2)*prmter[2];
        pattern4[1:16,3:4] = ones(16,2)*prmter[3];
        pattern4[1:16,5:6] = ones(16,2)*prmter[4];
        pattern4[1:16,7:8] = ones(16,2)*prmter[3];
        pattern4[1:16,9:10] = ones(16,2)*prmter[2];
        pattern4[1:16,11:12] = ones(16,2)*prmter[4];
        pattern4[1:16,13:14] = ones(16,2)*prmter[5];
        pattern4[1:16,15:16] = ones(16,2)*prmter[2];
        
          beta1[10:25,10:25] = pattern4;
          beta1[10:25,42:57] = pattern3;
          beta1[42:57,10:25] = pattern2;
          beta1[42:57,42:57] = pattern1;
          
          beta2[10:25,10:25] = pattern1;
          beta2[10:25,42:57] = pattern2;
          beta2[42:57,10:25] = pattern3;
          beta2[42:57,42:57] = pattern4;
          
          beta3[10:25,10:25] = pattern3;
          beta3[10:25,42:57] = pattern4;
          beta3[42:57,10:25] = pattern1;
          beta3[42:57,42:57 ]= pattern2;
}
if (pattern == 4)
{
  
  pattern1 = ones(12,12)*prmter[2];
  pattern1[1:12,3:4] = ones(12,2)*prmter[5];
  pattern1[1:12,5:6] = ones(12,2)*prmter[4];
  pattern1[1:12,7:8] = ones(12,2)*prmter[3];
  pattern1[1:12,9:1] = ones(12,2)*prmter[4];
  pattern1[1:12,11:12] = ones(12,2)*prmter[3];
  
 
    pattern2 = ones(12,12)*prmter[2];
    pattern2[1:12,7:12] = ones(12,6)*prmter[3];
    
   
      pattern3 = ones(12,12)*prmter[2];
      pattern3[1:12,5:8] = ones(12,4)*prmter[4];
      pattern3[1:12,9:12] = ones(12,4)*prmter[3];
      
    
        pattern4 = ones(12,12)*prmter[2];
        pattern4[1:12,4:6]= ones(12,3)*prmter[5];
        pattern4[1:12,7:9] = ones(12,3)*prmter[4];
        pattern4[1:12,10:12] = ones(12,3)*prmter[3];
        
        
          beta1(11:22,11:22) = t(pattern4);
          beta1(11:22,43:54) = t(pattern3);
          beta1(43:54,11:22) = t(pattern2);
          beta1(43:54,43:54) = t(pattern1);
          
          beta2(11:22,11:22) = t(pattern1);
          beta2(11:22,43:54) = t(pattern2);
          beta2(43:54,11:22) = t(pattern3);
          beta2(43:54,43:54) = t(pattern4);
          
          beta3(11:22,11:22) = t(pattern3);
          beta3(11:22,43:54) = t(pattern4);
          beta3(43:54,11:22) = t(pattern1);
          beta3(43:54,43:54) = t(pattern2);
}
if(pattern == 5)
{
  pattern1 = ones(18,18)*prmter[1];
  pattern1[1:5,1:18] = ones(5,18)*prmter[2];
  pattern1[6:10,1:5] = ones(5,5)*prmter[2];
  pattern1[14:18,1:18] = ones(5,18)*prmter[3];
  pattern1[9:13,14:18] = ones(5,5)*prmter[3];
  
  
    pattern2 = ones(18,18)*prmter[1];
    pattern2[1:4,1:18] = ones(4,18)*prmter[3];
    pattern2[5:8,1:4]= ones(4,4)*prmter[3];
    pattern2[15:18,1:18] = ones(4,18)*prmter[4];
    pattern2[11:14,15:18] = ones(4,4)*prmter[4];
    
  
      pattern3 = ones(18,18)*prmter[1];
      pattern3[1:3,1:18] = ones(3,18)*prmter[4];
      pattern3[4:6,1:3]= ones(3,3)*prmter[4];
      pattern3[16:18,1:18] = ones(3,18)*prmter[5];
      pattern3[13:15,16:18] = ones(3,3)*prmter[5];
      
        pattern4 = ones(18,18)*prmter[1];
        pattern4[1:2,1:18] = ones(2,18)*prmter[5];
        pattern4[3:4,1:2] = ones(2,2)*prmter[5];
        pattern4[17:18,1:18] = ones(2,18)*prmter[2];
        pattern4[15:16,17:18] = ones(2,2)*prmter[2];
        
        
          beta1[9:26,9:26] = pattern1;
          beta1[9:26,41:58] = pattern3;
          beta1[41:58,9:26] = pattern2;
          beta1[41:58,41:58] = pattern4;
          
          beta2[9:26,9:26] = pattern1;
          beta2[9:26,41:58] = pattern2;
          beta2[41:58,9:26] = pattern3;
          beta2[41:58,41:58] = pattern4;
          
          beta3[9:26,9:26] = pattern3;
          beta3[9:26,41:58] = pattern4;
          beta3[41:58,9:26] = pattern1;
          beta3[41:58,41:58] = pattern2;
}
if(pattern == 6)
{
    
    pattern1 = ones(18,18)*prmter[1];
    pattern1[1:5,1:18] = ones(5,18)*prmter[2];
    pattern1[6:10,1:5] = ones(5,5)*prmter[2];
    pattern1[14:18,1:18] = ones(5,18)*prmter[3];
    pattern1[9:13,14:18] = ones(5,5)*prmter[3];
    
  
      pattern2 = ones(18,18)*prmter[1];
      pattern2[1:5,1:18]= ones(5,18)*prmter[4];
      pattern2[6:10,1:5] = ones(5,5)*prmter[4];
      pattern2[14:18,1:18] = ones(5,18)*prmter[5];
      pattern2[9:13,14:18] = ones(5,5)*prmter[5];
      
    
        pattern3 = ones(18,18)*prmter[1];
        pattern3[1:5,1:18] = ones(5,18)*prmter[5];
        pattern3[6:10,1:5] = ones(5,5)*prmter[5];
        pattern3[14:18,1:18] = ones(5,18)*prmter[3];
        pattern3[9:13,14:18] = ones(5,5)*prmter[3];
        
       
          pattern4 = ones(18,18)*prmter[1];
          pattern4[1:5,1:18] = ones(5,18)*prmter[2];
          pattern4[6:10,1:5] = ones(5,5)*prmter[2];
          pattern4[14:18,1:18] = ones(5,18)*prmter[4];
          pattern4[9:13,14:18] = ones(5,5)*prmter[4];
          
          
          
            beta1[9:26,9:26] = pattern1;
            beta1[9:26,41:58] = t(pattern3);
            beta1[41:58,9:26] = pattern2;
            beta1[41:58,41:58] = t(pattern4);
            
            beta2[9:26,9:26] = t(pattern1);
            beta2[9:26,41:58] = pattern2;
            beta2[41:58,9:26] = t(pattern3);
            beta2[41:58,41:58] = t(pattern4);
            
            beta3[9:26,9:26] = t(pattern3);
            beta3[9:26,41:58] = t(pattern4);
            beta3[41:58,9:26] = pattern1;
            beta3[41:58,41:58] = pattern2;
  }
if (pattern == 7) 
{
  pattern1 = ones(18,18)*prmter[1];
  pattern1[1:5,1:18] = ones(5,18)*prmter[2];
  pattern1[6:10,1:5] = ones(5,5)*prmter[2];
  pattern1[14:18,1:18] = ones(5,18)*prmter[3];
  pattern1[9:13,14:18] = ones(5,5)*prmter[3];
  
  
    pattern2 = ones(18,18)*prmter[1];
    pattern2[1:5,1:18] = ones(5,18)*prmter[4];
    pattern2[6:10,1:5] = ones(5,5)*prmter[4];
    pattern2[14:18,1:18] = ones(5,18)*prmter[5];
    pattern2[9:13,14:18] = ones(5,5)*prmter[5];
    
   
      pattern3 = ones(18,18)*prmter[1];
      pattern3[1:5,1:18] = ones(5,18)*prmter[5];
      pattern3[6:10,1:5] = ones(5,5)*prmter[5];
      pattern3[14:18,1:18] = ones(5,18)*prmter[3];
      pattern3[9:13,14:18] = ones(5,5)*prmter[3];
      
     
        pattern4 = ones(18,18)*prmter[1];
        pattern4[1:5,1:18] = ones(5,18)*prmter[2];
        pattern4[6:10,1:5] = ones(5,5)*prmter[2];
        pattern4[14:18,1:18] = ones(5,18)*prmter[4];
        pattern4[9:13,14:18] = ones(5,5)*prmter[4];
        
        
     
          beta1[9:26,9:26] = pattern1;
          beta1[9:26,41:58] = t(pattern3);
          beta1[41:58,9:26] = pattern2;
          beta1[41:58,41:58] = t(pattern4);
          beta1 = beta1*0;
          
          beta2[9:26,9:26] = t(pattern1);
          beta2[9:26,41:58] = pattern2;
          beta2[41:58,9:26] = t(pattern3);
          beta2[41:58,41:58] = t(pattern4);
          
          beta3[9:26,9:26] = t(pattern3);
          beta3[9:26,41:58] = t(pattern4);
          beta3[41:58,9:26] = pattern1;
          beta3[41:58,41:58] = pattern2;
}
  beta = array(rep(0,64*64*3),c(64,64,3));
  beta[,,1] = beta1;
  beta[,,2] = beta2;
  beta[,,3] = beta3;
  
  roiIdsandBeta=list(RoiId1=0,RoiId2=0,RoiId3=0,RoiId4=0,RoiId5=0,Beta=0);
  for(ii in 1:5)
  {
  tempId = c();
  if (pattern == 7)  
  {
  for (jj in 2:3)
  {
  tempBeta = beta[,,jj];
  tempBeta = c(tempBeta);
  tempId1 = which(tempBeta==prmter[ii]);
  tempId = cbind(tempId,tempId1);
  }
  }
  else
  {
  for(jj in 1:3)
  {
  tempBeta = beta[,,jj];
  tempBeta = c(tempBeta);
  tempId1 = which(tempBeta==prmter[ii]);
  tempId = cbind(tempId,tempId1);
  }
  }       
  roiIdsandBeta[[ii]] = tempId;
  
}
  roiIdsandBeta[[6]] = beta;
  roiIdsandBeta ;
}