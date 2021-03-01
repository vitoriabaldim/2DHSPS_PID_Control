'''
    -> Proxima Temperatura      OK
    -> Controle Proporcional    OK
    -> Mudança de SetPoint      OK
    -> Controle Integral        OK
    -> Condições de Parada      OK
    -> Controle Derivativo      OK      
'''

#Bibliotecas
import numpy as np
import matplotlib.pyplot as plt

# --- FUNÇÃO PRINCIPAL ---
def main():
  P_total=40

  T_amb = 20.0 + 273.15
  T_alvo= 150.0 + 273.15

  Pot=0; T=T_amb

  dt =0.1; eps=1.7
  N_min=1200.0; N_max=18000.0

  kp=1.0
  ki=0.0
  kd=0.0

  VetorTempo=[]; VetorTemperatura=[T]; VetorSetPoint=[]
  VetorErros=[]; VetorErroArea = []; VetorIntegral = []
  VetorVarErro=[]; VetorPID = []; VetorPotencia = []

  cont=Simulacao_PID(VetorTempo,VetorTemperatura,VetorErros,VetorErroArea,VetorIntegral,VetorVarErro,VetorPID,VetorPotencia,VetorSetPoint,P_total,T_amb,T_alvo,T,Pot,eps,N_min,N_max,dt,kp,ki,kd)
  print(cont, "\nDuração da Simulação =", (cont)*dt, "s")

#  ExibiVetores(VetorTempo,VetorTemperatura,VetorErros,VetorErroArea,VetorIntegral,VetorVarErro,VetorPID,VetorPotencia,VetorSetPoint)

  Graficos(VetorTempo,VetorTemperatura,VetorErros,VetorErroArea,VetorIntegral,VetorVarErro,VetorPID,VetorPotencia,VetorSetPoint)

# --- FUNÇÃO QUE EXECUTA A SIMULAÇÃO DO PID ---
def Simulacao_PID(VetorTempo,VetorTemperatura,VetorErros,VetorErroArea,VetorIntegral,VetorVarErro,VetorPID,VetorPotencia,Referencia,P_total,T_amb,SetPoint,T,Pot,eps,N_min,N_max,dt,kp,ki,kd):
  i=0

  erro=Error(SetPoint, T)
  VetorErros.append(erro)

  var_erro=erro

  while( (((np.abs(erro)>=eps) or (np.abs(var_erro))>=eps) or (i<=N_min)) and (i<=N_max) ):                                               # AND: N funciona como um maximo de ciclos; OR: N funciona como um minimo de ciclos
    VetorTempo.append(i*dt)
    Referencia.append(SetPoint-273.15)

    #Medida da Temperatura e Calculo do Erro
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

    #Potencia de Saida
    Pot=PID(VetorPID,P_total,kp,ki,kd,erro,IntegralErro,var_erro)
    VetorPotencia.append(Pot)

    i+=1

    '''if(i==N_min/2):
      SetPoint=80+273.15'''

  return i-1

#Função que Calcula a Proxima Temperatura
def ProximaTemperatura(T,Pot,dt,T_amb):
  h=0.05; n=1.0; c=500.0; m_g=8.0
  C=(m_g/1000)*c
  T=T+((Pot/C)*(dt))-(h*dt*((T-T_amb)**n))
  return T

def Error(SetPoint, T_medido):
  erro=SetPoint-T_medido
  return erro

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

def ExibiVetores(VetorTempo,VetorTemperatura,VetorErros,VetorErroArea,VetorIntegral,VetorVarErro,VetorPID,VetorPotencia,VetorSetPoint):
  print("Tempo:",VetorTempo, '\n  -Quantidade de elementos:', len(VetorTempo), "; Ultimo elemento:",VetorTempo[-1])
  print("Temperatura:",VetorTemperatura, '\n  -Quantidade de elementos:', len(VetorTemperatura),"; Ultimo elemento:",VetorTemperatura[-1])
  print("Erros:",VetorErros, '\n  -Quantidade de elementos:', len(VetorErros),"; Ultimo elemento:",VetorErros[-1])
  print("Erro Area:",VetorErroArea, '\n  -Quantidade de elementos:', len(VetorErroArea),"; Ultimo elemento:",VetorErroArea[-1])
  print("Integral do Erro:",VetorIntegral, '\n  -Quantidade de elementos:', len(VetorIntegral),"; Ultimo elemento:",VetorIntegral[-1])
  print("Derivada do Erro:",VetorVarErro, '\n  -Quantidade de elementos:', len(VetorVarErro),"; Ultimo elemento:",VetorVarErro[-1])
  print("Potencia:",VetorPotencia, '\n  -Quantidade de elementos:', len(VetorPotencia),"; Ultimo elemento:",VetorPotencia[-1])

def Graficos(Tempo,Temperatura,Erro,ErroArea,Integral,Derivada,PID,Potencia,Referencia):
  del Temperatura[0]
  del Erro[0]
  
  TemperaturaC=[]
  ReferenciaZero=[]
  for i in range(len(Temperatura)):
    TemperaturaC.append(Temperatura[i]-273.15)
    ReferenciaZero.append(0)

  for i in range(len(Derivada)):
    if Derivada[i]>1500:
      Derivada[i]=0

  TemperaturaC1=TemperaturaC[900:1000]
  Referencia1=Referencia[900:1000]
  Tempo1=Tempo[900:1000]

#  plt.plot(Tempo,ReferenciaZero)
  plt.plot(Tempo1,Referencia1)
  plt.plot(Tempo1,TemperaturaC1)
#  plt.plot(Tempo,Erro)
#  plt.plot(Tempo,Integral)
#  plt.plot(Tempo,Derivada)
#  plt.plot(Tempo,PID)
#  plt.plot(Tempo,Potencia)
#  plt.plot(TemperaturaC,Potencia)

  desvio_padrao=np.std(TemperaturaC1)
  print(desvio_padrao)

  temp_media= np.mean(TemperaturaC1)
  print(temp_media)


# --- CHAMADA DA FUNÇÃO PRINCIPAL ---
main()
