# -*- coding:utf-8 -*-
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import PIL.Image as Image
import numpy as np
import math

omega_e = 7.2921151467e-5 # 地球自转速度
PI = 3.1415926535
GM = 3.986005e14

class Ephemeris(): # 星历类
    def __init__(self): # 各类参数
        self.AlreadyLoadData = 0 # 记录是否已经从数据文件中读取参数
        self.rA = 0
        self.C_rs = 0
        self.delta_n = 0
        self.M0 = 0
        self.C_us = 0
        self.e = 0
        self.C_uc = 0
        self.sqrt_a = 0
        self.t_oe = 0
        self.C_ic = 0
        self.omega = 0
        self.C_is = 0
        self.I0 = 0
        self.C_rc = 0
        self.w = 0
        self.omega_dot = 0
        self.I_dot = 0
        self.T_GD = 0
        self.t_oc = 0
        self.observt = 0
    
    def LoadEphemeris(self,filepath): # 从数据文件中读取参数
        if self.AlreadyLoadData: # 已经读取过的话则不再读取
            pass
        else:
            with open(filepath,"r") as fpread:
                data = fpread.readline().split(' ') # 数据在文件中以空格隔开
                self.C_rs = float(data[0])
                self.delta_n = float(data[1])
                self.M0 = float(data[2])
                self.C_uc = float(data[3])
                self.e = float(data[4])
                self.C_us = float(data[5])
                self.sqrt_a = float(data[6])
                self.t_oe = float(data[7])
                self.C_ic = float(data[8])
                self.omega = float(data[9])
                self.C_is = float(data[10])
                self.I0 = float(data[11])
                self.C_rc = float(data[12])
                self.w = float(data[13])
                self.omega_dot = float(data[14])
                self.I_dot = float(data[15])
                self.T_GD = float(data[16])
                self.t_oc = float(data[17])
            self.observt = self.t_oe # 设定观察时间从toe开始
            self.rA = self.sqrt_a * self.sqrt_a # 计算轨道半径
            self.AlreadyLoadData = 1 # 数据已读取

class Orbit_and_Satellite():
    def __init__(self): # 卫星和轨道参数
        self.satex = 5140678.179333044 # 缺省卫星位置
        self.satey = 26089623.83380793 # 缺省卫星位置
        self.I = PI*(55)/180 # 缺省轨道倾角
        self.L = PI*60/180 # 缺省轨道赤径
        self.e = 0.006164086284 # 缺省轨道偏心率
        self.rA = 5153.69580*5153.69580 # 缺省轨道半径
    
    def CalculateParameter(self, ephemeris): # 计算轨道参数和卫星位置
        if not ephemeris.AlreadyLoadData: # 如果没有读取过数据则无需计算，采用默认数据
            pass
        else:
            # 按照课程中讲授的顺序计算
            n0 = math.sqrt(GM/math.pow(ephemeris.rA,3))
            n = n0 + ephemeris.delta_n
            dt = ephemeris.observt - ephemeris.t_oe
            M = ephemeris.M0 + n * dt
            E = M
            E0 = 0

            while abs(E - E0) > 1e-6:
                E0 = E
                E = E0 - (E0 - M - ephemeris.e * math.sin(E0)) / (1 - ephemeris.e * math.cos(E0))

            f = math.atan2(math.sqrt(1 - math.pow(ephemeris.e, 2)) * math.sin(E), (math.cos(E) - ephemeris.e))
            u = ephemeris.w + f
            rt = ephemeris.rA * (1 - ephemeris.e * math.cos(E))
            I_t = ephemeris.I0 + ephemeris.I_dot * dt
            U = u + ephemeris.C_uc * math.cos(2 * u) + ephemeris.C_us * math.sin(2 * u)
            R = rt + ephemeris.C_rc * math.cos(2 * u) + ephemeris.C_rs * math.sin(2 * u)
            I = I_t + ephemeris.C_ic * math.cos(2 * u) + ephemeris.C_is * math.sin(2 * u)
            L = ephemeris.omega + ephemeris.omega_dot * dt

            # 存储绘图所需参数
            self.satex = R * math.cos(U)
            self.satey = R * math.sin(U)
            self.I = I
            self.L = L
            self.rA = ephemeris.rA
            self.e = ephemeris.e

    def drawSatellite(self): # 在卫星轨道平面内绘制卫星
        glPushMatrix() # 记录当前世界坐标系信息
        glTranslatef(self.satex/100000,self.satey/100000,0) # 将世界坐标系中心移到卫星所在位置
        # 在现在的世界坐标系中心画球，标识卫星
        glColor3f(0, 1, 1) # 设定卫星颜色
        sate = gluNewQuadric() # 初始化二次曲面变量
        gluSphere(sate,6,10,10) # 画圆，半径为6，后面的数字越大代表对圆的近似程度越大
        glColor3f(1, 1, 1) # 清空颜色，防止颜色对后面造成影响
        # 恢复世界坐标系
        glPopMatrix()
    
    def drawOrbit(self): # 在卫星轨道平面内绘制轨道
        glEnable(GL_BLEND)
        glEnable(GL_LINE_SMOOTH)
        slid = 100 # 轨道的近似多边形边数
        glLineWidth(2) # 设定线宽 
        glColor3f(0.2, 0.5, 1) # 设定线颜色
        glBegin(GL_LINE_LOOP) # 开始绘制
        #计算轨道参数
        ra = self.rA
        rc = ra*self.e
        rb = math.sqrt(ra*ra-rc*rc)
        theta = 0
        # 以多边形近似椭圆
        for i in range(slid):
            theta += 2*PI/slid
            x = ra*math.cos(theta)-rc
            y = rb*math.sin(theta)
            glVertex3f(x/100000,y/100000,0)
        glColor3f(1, 1, 1) # 清除颜色
        glEnd() # 结束绘制

class SatelliteSystem():
    def __init__(self):
        self.ephemeris = Ephemeris() # 星历
        self.orbit_and_satellite = Orbit_and_Satellite() # 轨道和卫星
        self.AlreadyLoadData = 0 # 记录是否已经读取数据

    def run(self):
        self.InitOpenGL() # 绘图初始化
        self.LoadTexture() # 加载纹理
        #调用函数绘制图像
        glutDisplayFunc(self.drawFunc) # 最开始默认执行一次drawFunc作为初始界面
        glutIdleFunc(self.drawFunc) # 设置 默认界面绘制函数 为 drawFunc
        glutKeyboardFunc(self.keyboard) # 设置 默认键盘响应函数 为 keyboard
        #主循环
        glutMainLoop() # 开始主循环

    def LoadTexture(self): # 加载纹理
        img = Image.open("earth.jpg").transpose(Image.FLIP_TOP_BOTTOM) # 打开图片（绝对不能是png
        img_data = np.array(list(img.getdata()),np.uint8) # 图片转成np数据
        texture = glGenTextures(1) #声明一个纹理标识
        glBindTexture(GL_TEXTURE_2D, texture) #绑定纹理标识
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) #设置纹理过滤模式
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) #设置纹理过滤模式
        glTexImage2D(GL_TEXTURE_2D,0,3,img.size[0],img.size[1],0,GL_RGB,GL_UNSIGNED_BYTE,img_data) #将图片绑定到纹理上     
        glEnable(GL_TEXTURE_2D) #启用纹理

    def InitOpenGL(self): # 显示界面初始化
        #使用glut初始化OpenGL
        glutInit()
        #显示模式:GLUT_SINGLE无缓冲直接显示|GLUT_RGBA采用RGB(A非alpha)
        glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH)#一定要添加深度检测缓冲区！！！也就是第三个参数
        #窗口位置及大小-生成
        glutInitWindowPosition(200,250)
        glutInitWindowSize(800,800)
        glutCreateWindow("卫星导航编程作业 by 小陈 using Python") # 窗体名字
        glViewport(0, 0, 1, 1)#设置显示范围
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45,1,10,10000)#设置视场
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glTranslatef(0, 0, -700)#平移
        #旋转坐标轴以便显示
        glRotatef(10, 1, 0 , 0)
        glRotatef(-30, 0, 1, 0)
        #将坐标轴旋转到正确位置
        glRotatef(-90, 1,0,0)
        glRotatef(-90, 0,0,1)

    def keyboard(self,key,x,y): # 键盘响应函数
        if key == b'r' and self.AlreadyLoadData == 0: # 如果按下r且未读取过数据则执行
            self.ephemeris.LoadEphemeris("data.txt")
            self.AlreadyLoadData = 1
    
    def drawFunc(self): # 绘制地球、赤道、地球坐标系、卫星和轨道
        self.drawInit() # 初始化
        glPushMatrix() # 记录当前世界坐标系
        # 将当前世界坐标系旋转为地球坐标系
        glRotatef(0.5*(self.ephemeris.observt-self.ephemeris.t_oe)/120,0,0,1) 
        self.drawAxis() # 绘制坐标轴
        self.drawEquator() # 绘制赤道
        self.drawEarth() # 绘制地球
        glPopMatrix() # 恢复世界坐标系
        self.drawOrbit_and_Satellite() # 绘制卫星和轨道
        glFlush() # 刷新显示
        self.ephemeris.observt += 120 # 观测时间点前进2min

    def drawInit(self):# 绘图初始化
        # 清除之前画面
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT) 
        # 启用深度检测
        glDepthMask(GL_TRUE)
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
        glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST)
        # 去除二次曲面的背面显示
        glEnable(GL_CULL_FACE)
        glShadeModel(GL_SMOOTH)
    
    def drawAxis(self): # 绘制坐标轴
        glEnable(GL_BLEND)
        glEnable(GL_LINE_SMOOTH)
        glLineWidth(2) # 设置线宽
        glBegin(GL_LINES) # 开始绘制
        # Z轴
        glColor3f(0, 0, 1)
        glVertex3f(0, 0, 0)
        glVertex3f(0, 0, 300)
        # X轴
        glColor3f(1, 0, 0)
        glVertex3f(0, 0, 0)
        glVertex3f(300, 0, 0)
        # Y轴
        glColor3f(0, 1, 0)
        glVertex3f(0, 0, 0)
        glVertex3f(0, 300, 0)
        glColor3f(1, 1, 1) # 清除颜色信息
        glEnd() # 结束绘制
        
    def drawEquator(self): # 绘制赤道
        glEnable(GL_BLEND)
        glEnable(GL_LINE_SMOOTH)
        slid = 100 # 以多边形近似椭圆
        glLineWidth(4) # 设定线宽
        glColor3f(1, 1, 0) # 设定颜色
        glBegin(GL_LINE_LOOP) # 开始绘制
        rd = 63.71 + 0 # 设定半径
        for i in range(slid):
            glVertex3f(rd *math.cos(2*PI / slid*i), rd*math.sin(2*PI / slid*i),0)
        glColor3f(1, 1, 1) # 清除颜色
        glEnd() # 结束绘制

    def drawEarth(self):
        glPushMatrix() # 记录当前世界坐标系
        glRotatef(90, 0, 0, 1) # 旋转世界坐标系以使纹理对齐格林威治本初子午线
        earth_model = gluNewQuadric() # 初始化二次曲面对象
        gluQuadricDrawStyle(earth_model, GL_FILL)
        gluQuadricNormals(earth_model, GLU_SMOOTH)
        gluQuadricTexture(earth_model, GL_TRUE) # 启用纹理
        gluSphere(earth_model,63.71,100,100) # 绘制地球
        glPopMatrix() # 恢复世界坐标系

    def drawOrbit_and_Satellite(self):
        self.orbit_and_satellite.CalculateParameter(self.ephemeris) # 计算参数
        glPushMatrix() # 记录当前世界坐标系
        # 旋转世界坐标系，使XOY平面与轨道平面重合
        glRotatef(180*self.orbit_and_satellite.L/PI,0,0,1)
        glRotatef(180*self.orbit_and_satellite.I/PI,1,0,0)
        self.orbit_and_satellite.drawOrbit() # 绘制轨道
        self.orbit_and_satellite.drawSatellite() # 绘制卫星
        glPopMatrix() # 恢复世界坐标系

if __name__ == '__main__':
    SS = SatelliteSystem()
    SS.run()