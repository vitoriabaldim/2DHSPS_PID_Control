#Bibliotecas
import numpy as np
import matplotlib.pyplot as plt
 
# --- FUNÇÃO PRINCIPAL ---
def main():
  #Parâmetros do Sistema e da Simulação
  P_total=40 
  T_amb = 20.0 + 273.15
  T_alvo= 40.0 + 273.15
  dt=0.01; eps=1.7
  N_min=9000.0; N_max=18000.0
 
  #Inicialização de Algumas Variáveis
  Pot=0; T=T_amb
 
  #Parâmetros do PID
  kp=1.0
  ki=0.5
  kd=0.01
 
  #Criação dos Vetores
  VetorTempo=[]; VetorTemperatura=[T]; VetorSetPoint=[]
  VetorErros=[]; VetorErroArea = []; VetorIntegral = []
  VetorVarErro=[]; VetorPID = []; VetorPotencia = []
 
  #Chamada da Função que executa a Simulação 
  Simulacao_PID(VetorTempo, VetorTemperatura, VetorErros, VetorErroArea, VetorIntegral, VetorVarErro, VetorPID, VetorPotencia, VetorSetPoint,P_total,T_amb,T_alvo,T,Pot,eps,N_min,N_max,dt,kp,ki,kd)
 
  #Chamada da Função que plota os gráficos
  Graficos(VetorTempo,VetorTemperatura,VetorSetPoint)
 
# --- FUNÇÃO QUE EXECUTA A SIMULAÇÃO DO PID ---
def Simulacao_PID(VetorTempo,VetorTemperatura,VetorErros, VetorErroArea,VetorIntegral,VetorVarErro,VetorPID,VetorPotencia,
Referencia,P_total,T_amb,SetPoint,T,Pot,eps,N_min,N_max,dt,kp,ki,kd):
  #inicialização de algumas variáveis para que se possa entrar no loop do while
  i=0
  erro=Error(SetPoint, T)
  VetorErros.append(erro)
  var_erro=erro
 
  #Loop: Simulação só para se: o erro se torna relativamente constante (variação < eps) e atingir um número mínimo de passos; ou então se atinge um número máximo de passos
  while( ((np.abs(var_erro)>=eps) or (i<=N_min)) and (i<=N_max) ):
    VetorTempo.append(i*dt)
    Referencia.append(SetPoint-273.15)
 
    #Medida da Temperatura e Cálculo do Erro
    T=ProximaTemperatura(T,Pot,dt,T_amb)
    VetorTemperatura.append(T)
    erro=Error(SetPoint, T)
    VetorErros.append(erro)
 
    #Integral do Erro
    erro_area=dt*((VetorErros[i]+VetorErros[i+1])/2)
    VetorErroArea.append(erro_area)
    IntegralErro=np.sum(VetorErroArea)
    VetorIntegral.append(IntegralErro)
 
    #Derivada do Erro
    var_erro=(VetorErros[i+1]-VetorErros[i])/dt
    VetorVarErro.append(var_erro)
 
    #Potência de Saída
    Pot=PID(VetorPID,P_total,kp,ki,kd,erro,IntegralErro,var_erro)
    VetorPotencia.append(Pot)
 
    #Incremento do contador i
    i+=1
 
    #Troca do SetPoint na metade da Simulação
    if(i==N_min/2):
      SetPoint=80+273.15
 
 
#Função que Calcula a Próxima Temperatura
def ProximaTemperatura(T,Pot,dt,T_amb):
  h=0.05; n=1.25; c=500.0; m_g=8.0
  C=(m_g/1000)*c
  T=T+((Pot/C)*(dt))-(h*dt*((T-T_amb)**n))
  return T
 
#Função que Calcula o Erro
def Error(SetPoint, T_medido):
  erro=SetPoint-T_medido
  return erro
 
#Função que calcula os parâmetros P, I, D e retorna a Potência de Saída
def PID(VetorPID,P_total,kp,ki,kd,erro,integral,derivada):
  P=kp*erro
  I=ki*integral
  D=kd*derivada
 
  PID=P+I+D
  VetorPID.append(PID)
 
  if(PID>1):
    Pot_Saida=P_total
  elif(PID<0):
    Pot_Saida=0
  else:
    Pot_Saida=P_total*PID
 
  return Pot_Saida
 
def Graficos(Tempo,Temperatura,Referencia):
  del Temperatura[0]
  
  TemperaturaC=[]
  for i in range(len(Temperatura)):
    TemperaturaC.append(Temperatura[i]-273.15)
 
  plt.plot(Tempo,Referencia)
  plt.plot(Tempo,TemperaturaC)
 
# --- CHAMADA DA FUNÇÃO PRINCIPAL ---
main()
