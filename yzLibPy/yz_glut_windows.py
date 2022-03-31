from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import numpy as np
import math

class GLUTWindowManager(object):
    def __init__(self):
        self.idle_func = []

    def AddIdleFunc(self, func):
        if callable(func):
            self.idle_func.append(func)

    def EnterMainLoop(self):
        glutIdleFunc(self.DefaultIdleFunc)
        glutMainLoop()

    def DefaultIdleFunc(self):
        for func in self.idle_func:
            if callable(func):
                func()


class BaseGLUTWindow(object):
    def __init__(self):
        self.win_id = 0
        self.win_x = 0
        self.win_y = 0
        self.win_width = 800
        self.win_height = 600

        self.old_x = 0
        self.old_y = 0
        self.mouse_state = [GLUT_UP, GLUT_UP, GLUT_UP]
        self.modifier = 0

    def RegisterOpenGLCallbacks(self):
        funcs = [func for func in dir(self) if callable(getattr(self, func))]
        if 'DisplayFunc' in funcs:
            glutDisplayFunc(getattr(self, 'DisplayFunc'))
        if 'ReshapeFunc' in funcs:
            glutReshapeFunc(getattr(self, 'ReshapeFunc'))
        if 'KeyboardFunc' in funcs:
            glutKeyboardFunc(getattr(self, 'KeyboardFunc'))
        if 'MouseFunc' in funcs:
            glutMouseFunc(getattr(self, 'MouseFunc'))
        if 'MotionFunc' in funcs:
            glutMotionFunc(getattr(self, 'MotionFunc'))
        if 'PassiveMotionFunc' in funcs:
            glutPassiveMotionFunc(getattr(self, 'PassiveMotionFunc'))
        if 'SpecialFunc' in funcs:
            glutSpecialFunc(getattr(self, 'SpecialFunc'))
        if 'IdleFunc' in funcs:
            glutIdleFunc(getattr(self, 'IdleFunc'))


class GLUTWindow3D(BaseGLUTWindow):
    def __init__(self):
        super().__init__()

        self.eye_x = 0.0
        self.eye_y = 0.0
        self.eye_z = 3.0
        self.spin_x = 0.0
        self.spin_y = 0.0
        self.fovy = 45.0
        self.z_near = 0.1
        self.z_far = 100.0
        self.use_arcball_flag = False   # currently not ready for work
        self.arcball_mat_col_major = np.identity(4).astype(dtype=float)

        self.back_ground_red = 1.0
        self.back_ground_green = 1.0
        self.back_ground_blue = 1.0
        self.back_ground_alpha = 1.0

    def CreateGLUTWindow(self):
        glutInit()
        glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA)
        glutInitWindowPosition(self.win_x, self.win_y)
        glutInitWindowSize(self.win_width, self.win_height)
        self.win_id = glutCreateWindow(b"3D Window")

        self.RegisterOpenGLCallbacks()

    def Draw(self):
        glColor3f(1, 0, 0)
        glutSolidTeapot(0.5)

    def DrawAppend(self):
        pass

    def SetLightMaterial(self):
        ambient = (0.0, 0.0, 0.0, 1.0)
        diffuse = (1.0, 1.0, 1.0, 1.0)
        specular = (1.0, 1.0, 1.0, 1.0)
        position = (10.0, 20.0, 5.0, 0.0)

        material_ambient = (0.1, 0.1, 0.1, 1.0)
        material_diffuse = (0.5, 0.5, 0.5, 1.0)
        material_specular = (0.0, 0.0, 0.0, 1.0)
        shineness = 0.0

        # setup light
        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular)

        # set light position, we set the light to be rotation free with mouse
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        glLightfv(GL_LIGHT0, GL_POSITION, position)
        glPopMatrix()

        # set material
        glEnable(GL_COLOR_MATERIAL)
        glMaterialfv(GL_FRONT, GL_AMBIENT, material_ambient)
        glMaterialfv(GL_FRONT, GL_DIFFUSE, material_diffuse)
        glMaterialfv(GL_FRONT, GL_SPECULAR, material_specular)
        glMaterialf(GL_FRONT, GL_SHININESS, shineness)

    def DisplayFunc(self):
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glShadeModel(GL_SMOOTH)

        glClearColor(self.back_ground_red, self.back_ground_green, self.back_ground_blue, self.back_ground_alpha)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # Setup Modelview Matrix
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        gluLookAt(self.eye_x, self.eye_y, self.eye_z, self.eye_x, self.eye_y, self.eye_z - 1, 0, 1, 0)

        if self.use_arcball_flag:
            glMultMatrixd(self.arcball_mat_col_major)
        else:   # use turn table
            glRotatef(self.spin_y, 1, 0, 0)
            glRotatef(self.spin_x, 0, 1, 0)

        self.SetLightMaterial()
        self.Draw()
        self.DrawAppend()

        glutSwapBuffers()

    def ReshapeFunc(self, w, h):
        self.win_width = w
        self.win_height = h

        # Setup Projection Matrix
        glViewport(0, 0, self.win_width, self.win_height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(self.fovy, float(self.win_width) / float(self.win_height), self.z_near, self.z_far)

        # set back to model view
        glMatrixMode(GL_MODELVIEW)

    def KeyboardFunc(self, key, x, y):
        if key == b'\x1b':      # ESC, ascii 27
            exit(0)

        print(key)

        glutPostRedisplay()

    def MouseFunc(self, button, state, x, y):
        if button < 0 or button > 2:
            return

        self.mouse_state[button] = state
        self.old_x = x
        self.old_y = y
        self.modifier = glutGetModifiers()

        glutPostRedisplay()

    def MotionFunc(self, x, y):
        if self.mouse_state[GLUT_LEFT_BUTTON] == GLUT_DOWN:
            if self.modifier & GLUT_ACTIVE_SHIFT:
                self.eye_x -= float(x - self.old_x) / 20
                self.eye_y += float(y - self.old_y) / 20
            else:
                if self.use_arcball_flag:
                    if x != self.old_x or y != self.old_y:
                        va = self.GetArcballVector(self.old_x, self.old_y)
                        vb = self.GetArcballVector(x, y)
                        angle = math.acos(min(1.0, np.inner(va, vb)))
                        axis_in_camera_coord = np.cross(va, vb)
                else:
                    self.spin_x += float(x - self.old_x)
                    self.spin_y += float(y - self.old_y)
        elif self.mouse_state[GLUT_MIDDLE_BUTTON] == GLUT_DOWN:
            self.eye_x -= float(x - self.old_x) / 20
            self.eye_y += float(y - self.old_y) / 20
        elif self.mouse_state[GLUT_RIGHT_BUTTON] == GLUT_DOWN:
            self.eye_z += float(y - self.old_y) / 10

        self.old_x = x
        self.old_y = y
        glutPostRedisplay()

    def IdleFunc(self):
        glutSetWindow(self.win_id)
        glutPostRedisplay()

    def GetArcballVector(self, x, y):
        px = x / float(self.win_width) * 2 - 1.0
        py = - (y / float(self.win_height) * 2 - 1.0)
        op_square = px * px + py * py
        if op_square <= 1.0:
            pz = math.sqrt(1.0 - op_square)
        else:
            norm = math.sqrt(op_square)
            px /= norm
            py /= norm
        return np.array([px, py, pz]).astype(dtype=float)

